"""
tests/physics/test_rk3.py
--------------------------
Tests specific to the SSP-RK3 time integration scheme (Shu & Osher 1988).

Four test classes corresponding to the four properties identified after
replacing the Hancock predictor-corrector with SSP-RK3:

  1. TestSSPProperty
       No new extrema should appear in smooth problems (field loop, OT).
       The SSP property guarantees this when the spatial operator is TVD.

  2. TestRK3EnergyConservation
       Total energy loss at t=0.5 on the OT vortex must be < 3%.
       RK3 reduces temporal truncation error relative to Hancock so this
       bound should be comfortably achievable.

  3. TestRK3StageCoefficients
       Mass must be conserved to machine precision for periodic BCs.
       Any error in the blend weights (3/4+1/4 or 1/3+2/3) introduces a
       systematic mass source or sink that this regression catches.
"""

import numpy as np
import h5py
import pytest
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _energy_history(h5_path):
    """
    Return (times, E_mean) where E_mean is the domain-averaged total energy
    density (thermal + kinetic + magnetic) at each snapshot.

    Uses gamma from parameters/gamma in the HDF5.  The gamma_actual = gamma
    fix in main.f90 ensures this is correct for all problem types.
    """
    with h5py.File(h5_path, "r") as f:
        gamma = float(np.asarray(f["parameters/gamma"]).item())
        if gamma <= 1.0:
            pytest.skip(
                "gamma <= 1 in HDF5 (gamma_actual not set for this problem type). "
                "Add 'gamma_actual = gamma' after initialize_problem() in main.f90."
            )
        times = np.array(f["time/sim_time"])
        keys  = sorted(f["rho"].keys(), key=lambda k: int(k.replace("SNAPSHOT", "")))
        E_mean = []
        for k in keys:
            rho = np.array(f["rho"][k])
            P   = np.array(f["P"][k])      # thermal pressure (magnetic term stripped)
            Bx  = np.array(f["Bx"][k])
            By  = np.array(f["By"][k])
            Vx  = np.array(f["Vx"][k])
            Vy  = np.array(f["Vy"][k])
            e   = (P / (gamma - 1.0)
                   + 0.5 * rho * (Vx**2 + Vy**2)
                   + 0.5 * (Bx**2 + By**2))
            E_mean.append(float(e.mean()))
    return times, np.array(E_mean)


def _mass_history(h5_path):
    """Return the total (unnormalised) density sum at each snapshot."""
    with h5py.File(h5_path, "r") as f:
        keys = sorted(f["rho"].keys(), key=lambda k: int(k.replace("SNAPSHOT", "")))
        return [float(np.array(f["rho"][k]).sum()) for k in keys]


# ---------------------------------------------------------------------------
# 1. SSP property
# ---------------------------------------------------------------------------

class TestSSPProperty:
    """
    SSP-RK3 should be non-oscillatory when the spatial operator is TVD.
    The slope limiter in CHIMERA makes the spatial operator TVD, so no new
    extrema should appear in smooth flows.

      - Field loop: advection, density should be preserved
      - OT vortex:  compressive; some growth, no large spikes
    """

    def test_field_loop_no_density_spike(self, run_sim):
        """
        For uniform-density advection (rho=1 everywhere), max(rho) must
        not grow by more than 0.1%.
        """
        h5 = run_sim(problem_type=3, N=64, h5_name="rk3_ssp_fl.h5", timeout=120)

        with h5py.File(h5, "r") as f:
            keys   = sorted(f["rho"].keys(), key=lambda k: int(k.replace("SNAPSHOT", "")))
            rho_0  = float(np.array(f["rho"][keys[0]]).max())
            rho_hi = max(float(np.array(f["rho"][k]).max()) for k in keys)

        tol = 1e-3 * rho_0
        assert rho_hi <= rho_0 + tol, (
            f"Field loop: max(rho) grew {100*(rho_hi - rho_0)/rho_0:.4f}% "
            f"(from {rho_0:.6f} to {rho_hi:.6f}, tolerance 0.1%). "
            f"SSP-RK3 with TVD reconstruction must not create new density maxima."
        )

    def test_ot_vortex_pressure_floor(self, run_sim):
        """
        OT vortex: minimum thermal pressure must never touch P_floor = 1e-12
        in more than a handful of cells.

        SSP-RK3 with TVD reconstruction is dissipative -- it should not create
        low-pressure rarefaction spikes. If many cells sit at P_floor, the scheme
        is driving thermal energy negative and the floors are masking an
        instability rather than catching an isolated near-vacuum cell.

        Threshold: at most 0.1% of cells (4 cells at N=64) at P_floor at any
        snapshot. A clean SSP scheme should have zero such cells on OT.
        """
        h5 = run_sim(problem_type=1, N=64, h5_name="rk3_ssp_ot.h5", timeout=120)

        P_floor = 1e-12
        max_floor_fraction = 0.0

        with h5py.File(h5, "r") as f:
            keys = sorted(f["P"].keys(), key=lambda k: int(k.replace("SNAPSHOT", "")))
            N_cells = np.array(f["P"][keys[0]]).size
            for k in keys:
                P = np.array(f["P"][k])
                fraction = float(np.sum(P <= P_floor * 10)) / N_cells
                max_floor_fraction = max(max_floor_fraction, fraction)

        assert max_floor_fraction < 1e-3, (
            f"OT vortex: {max_floor_fraction*100:.3f}% of cells at or near P_floor "
            f"at some snapshot (threshold 0.1%). SSP-RK3 should not drive thermal "
            f"pressure to the floor across multiple cells. This indicates rarefaction "
            f"instability, not physical behaviour."
        )



# ---------------------------------------------------------------------------
# 2. Energy conservation
# ---------------------------------------------------------------------------

class TestRK3EnergyConservation:
    """
    Ideal MHD should conserve total energy (only numerical dissipation).
    SSP-RK3 (3rd-order in time) reduces temporal dissipation relative to
    Hancock (2nd-order), so the energy loss bound should be tighter.

    Threshold: relative change in Energy < 3% at t=0.5 on the OT vortex at N=64..
    """

    def test_ot_energy_loss_at_half_period(self, run_sim):
        h5 = run_sim(problem_type=1, N=64, h5_name="rk3_energy_ot.h5", timeout=120)

        times, E_mean = _energy_history(h5)

        idx = int(np.argmin(np.abs(times - 0.5)))
        if abs(times[idx] - 0.5) > 0.02:
            pytest.skip(f"No snapshot near t=0.5; closest is t={times[idx]:.4f}")

        dE = abs((E_mean[idx] - E_mean[0]) / E_mean[0])

        assert dE < 0.03, (
            f"OT vortex: |ΔE/E₀| = {dE*100:.3f}% at t≈{times[idx]:.3f} "
            f"(threshold 3%). SSP-RK3 should conserve energy better than Hancock "
        )

    def test_energy_monotone_early(self, run_sim):
        """
        Total energy shouldn't go up during the early smooth phase (t < 0.2).
        Any increase is unphysical and indicates the scheme is adding energy
        """
        h5 = run_sim(problem_type=1, N=64, h5_name="rk3_energy_ot.h5", timeout=120)

        times, E_mean = _energy_history(h5)

        early = E_mean[times <= 0.2]
        if len(early) < 2:
            pytest.skip("Not enough early snapshots (t <= 0.2) to check monotonicity.")

        max_increase = max(0.0, float(np.max(np.diff(early))))

        assert max_increase < 1e-6 * abs(E_mean[0]), (
            f"Total energy increased by {max_increase:.2e} during t <= 0.2. "
            f"A TVD scheme with SSP-RK3 shouldn't add energy. This may indicate "
            f"anti-diffusion from incorrect RK3 blend signs."
        )


# ---------------------------------------------------------------------------
# 3. Stage coefficient regression
# ---------------------------------------------------------------------------

class TestRK3StageCoefficients:
    """
    For periodic BCs, the divergence theorem guarantees mass conservation
    regardless of time integration order.

    Any error in the SSP-RK3 blend coefficients (3/4 + 1/4 .ne. 1, or
    1/3 + 2/3 .ne. 1) introduces a mass source or sink.
    """

    def test_mass_conserved_ot_vortex(self, run_sim):
        """OT vortex is fully periodic — mass must be conserved to 1e-10."""
        h5 = run_sim(problem_type=1, N=64, h5_name="rk3_mass_ot.h5", timeout=120)

        masses    = _mass_history(h5)
        mass_0    = masses[0]
        max_drift = max(abs(m - mass_0) / mass_0 for m in masses)

        assert max_drift < 1e-10, (
            f"OT vortex mass drift: {max_drift:.2e} (threshold 1e-10). "
            f"For periodic BCs, mass must be conserved to machine precision. "
            f"A larger drift indicates an error in the SSP-RK3 blend weights "
            f"in main.f90 (check 3/4+1/4=1 and 1/3+2/3=1 in stages 2 and 3)."
        )

    def test_mass_conserved_field_loop(self, run_sim):
        """Field loop is periodic; use as a second check."""
        h5 = run_sim(problem_type=3, N=64, h5_name="rk3_mass_fl.h5", timeout=120)

        masses    = _mass_history(h5)
        mass_0    = masses[0]
        max_drift = max(abs(m - mass_0) / mass_0 for m in masses)

        assert max_drift < 1e-10, (
            f"Field loop mass drift: {max_drift:.2e} (threshold 1e-10). "
            f"Same diagnostic as OT vortex on a second periodic problem."
        )

    def test_mass_conserved_kh_instability(self, run_sim):
        """
        KH instability is periodic; use as a a third check with new
        physics (shear flow, no initial B field).
        """
        h5 = run_sim(problem_type=2, N=64, h5_name="rk3_mass_kh.h5", timeout=120)

        masses    = _mass_history(h5)
        mass_0    = masses[0]
        max_drift = max(abs(m - mass_0) / mass_0 for m in masses)

        assert max_drift < 1e-10, (
            f"KH instability mass drift: {max_drift:.2e} (threshold 1e-10)."
        )
