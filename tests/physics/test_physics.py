"""
Physics fidelity tests: verify that CHIMERA produces physically correct
results on the standard periodic test problems.

Checks:
  - Mass is conserved to machine precision (periodic BCs, no sources)
  - Density and thermal pressure are positive everywhere
  - Kelvin-Helmholtz completes without crash
  - Field loop advection converges at second order in resolution
"""

import numpy as np
import h5py

N_FAST = 32
N_CONV = 64   # for convergence ratio: compare N_FAST vs N_CONV


def _rho_snapshots(h5_path):
    """Return all rho snapshots stacked as (n_snap, N, N)."""
    with h5py.File(h5_path, "r") as f:
        grp  = f["rho"]
        keys = sorted(grp.keys(),
                      key=lambda k: int(k.replace("SNAPSHOT", "")))
        return np.stack([np.array(grp[k]) for k in keys])


def _field_magnitudes(h5_path, snap_index):
    """Return |B| at the given snapshot index (0-based)."""
    with h5py.File(h5_path, "r") as f:
        grp  = f["Bx"]
        keys = sorted(grp.keys(),
                      key=lambda k: int(k.replace("SNAPSHOT", "")))
        key  = keys[snap_index]
        Bx = np.array(f["Bx"][key])
        By = np.array(f["By"][key])
    return np.sqrt(Bx**2 + By**2)


class TestMassConservation:
    """
    For periodic BCs with no sources, total mass integral rho*dV must be
    conserved to machine precision throughout the simulation.
    """

    def test_ot_vortex_mass_conserved(self, run_sim):
        rho = _rho_snapshots(run_sim(problem_type=1, N=N_FAST))
        mass_first = rho[0].sum()
        mass_last  = rho[-1].sum()
        rel_err = abs(mass_last - mass_first) / mass_first
        assert rel_err < 1e-10, \
            f"OT vortex mass not conserved: relative error = {rel_err:.2e}"

    def test_field_loop_mass_conserved(self, run_sim):
        rho = _rho_snapshots(run_sim(problem_type=3, N=N_FAST))
        mass_first = rho[0].sum()
        mass_last  = rho[-1].sum()
        rel_err = abs(mass_last - mass_first) / mass_first
        assert rel_err < 1e-10, \
            f"Field loop mass not conserved: relative error = {rel_err:.2e}"


class TestFieldPositivity:
    """Density and thermal pressure must be positive everywhere, always."""

    def test_density_positive_ot(self, run_sim):
        rho = _rho_snapshots(run_sim(problem_type=1, N=N_FAST))
        assert rho.min() > 0.0, \
            f"Non-positive density in OT vortex: min = {rho.min():.3e}"

    def test_pressure_positive_ot(self, run_sim):
        with h5py.File(run_sim(problem_type=1, N=N_FAST), "r") as f:
            grp  = f["P"]
            keys = sorted(grp.keys(),
                          key=lambda k: int(k.replace("SNAPSHOT", "")))
            P_all = np.stack([np.array(grp[k]) for k in keys])
        assert P_all.min() > 0.0, \
            f"Non-positive thermal pressure in OT vortex: min = {P_all.min():.3e}"

    def test_density_positive_field_loop(self, run_sim):
        rho = _rho_snapshots(run_sim(problem_type=3, N=N_FAST))
        assert rho.min() > 0.0, \
            f"Non-positive density in field loop: min = {rho.min():.3e}"


class TestProblemCompletion:
    """All standard problems should complete without crashing."""

    def test_kh_instability_completes(self, run_sim):
        with h5py.File(run_sim(problem_type=2, N=N_FAST), "r") as f:
            n = len(f["rho"].keys())
        assert n >= 5, f"KH: only {n} snapshots written - may have crashed"


class TestFieldLoopConvergence:
    """
    Second-order convergence of the field loop advection scheme.

    The loop (vx=-2, vy=-1) returns to its starting position exactly at
    tEnd=1.0. The L1 error in |B| between the final and initial snapshots
    should decrease by a factor >= 3 when N doubles.

    A ratio < 3.0 (rather than the theoretical 4.0) is used to allow for
    slope-limiter effects near the loop edge, which degrade local accuracy
    at coarse resolution.
    """

    def _fla_l1_error(self, run_sim, N):
        """Normalised L1 error of |B| between t=tEnd and t=0."""
        h5 = run_sim(problem_type=3, N=N, h5_name=f"fla_conv_N{N}.h5")
        B0 = _field_magnitudes(h5, snap_index=0)
        Bf = _field_magnitudes(h5, snap_index=-1)
        mean_B0 = B0[B0 > 1e-8].mean() if (B0 > 1e-8).any() else 1.0
        return np.mean(np.abs(Bf - B0)) / mean_B0

    def test_error_decreases_with_resolution(self, run_sim):
        err_32 = self._fla_l1_error(run_sim, N_FAST)
        err_64 = self._fla_l1_error(run_sim, N_CONV)
        assert err_64 < err_32, (
            f"Error did not decrease with resolution: "
            f"N={N_FAST}: {err_32:.3e}, N={N_CONV}: {err_64:.3e}"
        )

    def test_second_order_convergence_rate(self, run_sim):
        err_32 = self._fla_l1_error(run_sim, N_FAST)
        err_64 = self._fla_l1_error(run_sim, N_CONV)
        ratio  = err_32 / err_64
        assert ratio >= 3.0, (
            f"Convergence ratio = {ratio:.2f} (expected >= 3.0 for second order). "
            f"L1 errors: N={N_FAST}: {err_32:.3e}, N={N_CONV}: {err_64:.3e}"
        )
