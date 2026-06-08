"""
Stress tests: verify that the solver behaves correctly under challenging
initial conditions.

Three regimes are tested:

1. Standard GRF (random sampling, seed fixed for reproducibility)
   - Checks completion and plausible parameter values.

2. High Mach number  (target_M_s = 5.0, > default max 2.0)
   - Exercises the speed-cap code path in setup_GRF_fields.
   - Verifies numerical stability and field positivity despite high
     initial velocities.

3. Low plasma beta  (target_beta = 0.1, < default min 0.5)
   - Creates a magnetically dominated initial state where thermal
     pressure is a small residual of kinetic + magnetic energy.
   - Exercises the proportional pressure floor (Issue 2 fix) and
     the MOOD reconstruction in the presence of large B.
"""

import numpy as np
import h5py

GRF_N    = 64
GRF_SEED = 42

# Plausibility bounds for the standard GRF (very generous; catch NaN/Inf/sign errors)
M_S_BOUNDS  = (0.01, 20.0)
BETA_BOUNDS = (0.01, 100.0)

# Stress-test targets (both outside the default random sampling range)
HIGH_MACH_TARGET = 5.0   # >> default max M_s = 2.0
LOW_BETA_TARGET  = 0.1   # << default min beta = 0.5


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _sorted_keys(grp):
    return sorted(grp.keys(), key=lambda k: int(k.replace("SNAPSHOT", "")))


def _min_over_all_snaps(h5_path, field):
    """Return the global minimum of `field` across all snapshots."""
    with h5py.File(h5_path, "r") as f:
        return min(np.array(f[field][k]).min() for k in _sorted_keys(f[field]))


def _snap_count(h5_path):
    with h5py.File(h5_path, "r") as f:
        return len(f["rho"].keys())


def _params(h5_path):
    """Return stored parameters dict.

    mhd_write_h5.f90 writes scalars as rank-1 shape-(1,) datasets
    (``scalar_dims(1) = 1`` + ``h5screate_simple_f(1, ...)``), not 0-d
    scalars.  ``np.asarray(...).item()`` extracts a Python scalar
    correctly regardless of whether the dataset is shape () or shape (1,).
    """
    with h5py.File(h5_path, "r") as f:
        return {
            "gamma": np.asarray(f["parameters/gamma"]).item(),
            "M_s":   np.asarray(f["parameters/M_s"]).item(),
            "beta":  np.asarray(f["parameters/beta"]).item(),
        }


# ---------------------------------------------------------------------------
# 1. Standard GRF
# ---------------------------------------------------------------------------

class TestGRFCompletion:
    """GRF runs should complete without crashing at any valid seed."""

    def test_grf_completes_default_seed(self, run_sim):
        n = _snap_count(run_sim(problem_type=5, N=GRF_N, seed=GRF_SEED))
        assert n >= 5, f"GRF seed={GRF_SEED}: only {n} snapshots written"

    def test_grf_completes_different_seed(self, run_sim):
        """Check a second seed to ensure robustness to random variation."""
        n = _snap_count(run_sim(problem_type=5, N=GRF_N, seed=GRF_SEED + 1))
        assert n >= 5, f"GRF seed={GRF_SEED + 1}: only {n} snapshots written"


class TestGRFParameters:
    """Realised diagnostics must be positive, finite, and physically plausible."""

    def test_grf_parameters_plausible(self, run_sim):
        """
        All three stored parameters (gamma, M_s, beta) checked in a single
        run to avoid spawning redundant simulations.
        """
        h5 = run_sim(problem_type=5, N=GRF_N, seed=GRF_SEED)
        p  = _params(h5)

        assert np.isfinite(p["gamma"]) and p["gamma"] > 1.0, \
            f"gamma invalid: {p['gamma']}"

        assert np.isfinite(p["M_s"]) and p["M_s"] > 0.0, \
            f"M_s invalid: {p['M_s']}"
        assert M_S_BOUNDS[0] < p["M_s"] < M_S_BOUNDS[1], \
            f"M_s = {p['M_s']:.3f} outside plausible range {M_S_BOUNDS}"

        assert np.isfinite(p["beta"]) and p["beta"] > 0.0, \
            f"beta invalid: {p['beta']}"
        assert BETA_BOUNDS[0] < p["beta"] < BETA_BOUNDS[1], \
            f"beta = {p['beta']:.3f} outside plausible range {BETA_BOUNDS}"

    def test_grf_fields_positive_throughout(self, run_sim):
        """Density and thermal pressure must be positive in all snapshots."""
        h5      = run_sim(problem_type=5, N=GRF_N, seed=GRF_SEED)
        min_rho = _min_over_all_snaps(h5, "rho")
        min_p   = _min_over_all_snaps(h5, "P")
        assert min_rho > 0.0, f"Non-positive density in GRF: min = {min_rho:.3e}"
        assert min_p   > 0.0, f"Non-positive pressure in GRF: min = {min_p:.3e}"


# ---------------------------------------------------------------------------
# 2. High Mach number stress test
# ---------------------------------------------------------------------------

class TestHighMachStress:
    """
    Stress test at target_M_s = 5.0 (> default cap of 2.0).

    At this target, setup_GRF_fields first scales velocities to RMS ~ V_0
    (large), then the speed-cap at c_f_max = 3 clips them back down.
    The initial condition therefore probes the speed-cap code path and
    starts the simulation in a state with significant velocity gradients.

    Because the speed cap reduces |v| substantially, the realised M_s_actual
    stored in the HDF5 will be lower than 5.0 (typically 1.5–2.5 depending
    on gamma and beta).  The assertion therefore checks numerical stability
    (completion + positivity) rather than an exact M_s value.
    """

    def test_high_mach_stable(self, run_sim):
        """
        Solver must complete and keep rho, p > 0 when initialised with
        target_M_s = 5.0 (speed-cap path).
        """
        h5 = run_sim(
            problem_type=5, N=GRF_N, seed=GRF_SEED,
            h5_name="stress_high_mach.h5",
            extra_args=[HIGH_MACH_TARGET, -1.0],   # M_s override, beta random
        )
        n       = _snap_count(h5)
        min_rho = _min_over_all_snaps(h5, "rho")
        min_p   = _min_over_all_snaps(h5, "P")

        assert n >= 5, \
            f"High-Mach (target M_s={HIGH_MACH_TARGET}): only {n} snapshots — possible crash"
        assert min_rho > 0.0, \
            f"Non-positive density at target M_s={HIGH_MACH_TARGET}: min = {min_rho:.3e}"
        assert min_p > 0.0, \
            f"Non-positive thermal pressure at target M_s={HIGH_MACH_TARGET}: min = {min_p:.3e}"

    def test_high_mach_parameters_finite(self, run_sim):
        """Stored diagnostics must be finite even after the speed-cap rescaling."""
        h5 = run_sim(
            problem_type=5, N=GRF_N, seed=GRF_SEED,
            h5_name="stress_high_mach_params.h5",
            extra_args=[HIGH_MACH_TARGET, -1.0],
        )
        p = _params(h5)
        for name, val in p.items():
            assert np.isfinite(val) and val > 0.0, \
                f"Parameter '{name}' = {val} is not finite/positive after high-Mach init"


# ---------------------------------------------------------------------------
# 3. Low plasma beta stress test
# ---------------------------------------------------------------------------

class TestLowBetaStress:
    """
    Stress test at target_beta = 0.1 (< default minimum of 0.5).

    In the magnetically dominated regime (beta << 1) the magnetic
    pressure B^2/2 >> p_thermal, so thermal pressure is a small
    residual of large numbers — exactly the cancellation scenario that
    the Issue 2 proportional pressure floor was designed to handle.

    Unlike the high-Mach case, the speed cap does NOT significantly
    alter the magnetic field (B is only rescaled in step 5, not in the
    speed cap), so beta_actual stored in the file should closely match
    target_beta.
    """

    def test_low_beta_stable(self, run_sim):
        """
        Solver must complete and keep rho, p > 0 in the magnetically
        dominated regime (target_beta = 0.1).
        """
        h5 = run_sim(
            problem_type=5, N=GRF_N, seed=GRF_SEED,
            h5_name="stress_low_beta.h5",
            extra_args=[-1.0, LOW_BETA_TARGET],    # M_s random, beta override
        )
        n       = _snap_count(h5)
        min_rho = _min_over_all_snaps(h5, "rho")
        min_p   = _min_over_all_snaps(h5, "P")

        assert n >= 5, \
            f"Low-beta (target beta={LOW_BETA_TARGET}): only {n} snapshots — possible crash"
        assert min_rho > 0.0, \
            f"Non-positive density at target beta={LOW_BETA_TARGET}: min = {min_rho:.3e}"
        assert min_p > 0.0, \
            f"Non-positive thermal pressure at target beta={LOW_BETA_TARGET}: min = {min_p:.3e}"

    def test_low_beta_realized(self, run_sim):
        """
        The realised plasma beta stored in the HDF5 must be in the
        magnetically dominated regime (< 0.5).

        Derivation: beta_actual = 2*P_avg / B_rms^2 is computed before
        the total-pressure conversion in setup_GRF_fields, so
        P_avg ≈ P_mean = 1.0 and B_rms ≈ B_0 = sqrt(2/target_beta),
        giving beta_actual ≈ target_beta regardless of the random seed.
        """
        h5 = run_sim(
            problem_type=5, N=GRF_N, seed=GRF_SEED,
            h5_name="stress_low_beta_val.h5",
            extra_args=[-1.0, LOW_BETA_TARGET],
        )
        beta_actual = _params(h5)["beta"]
        assert beta_actual < 0.5, (
            f"Realised beta = {beta_actual:.4f} is not magnetically dominated "
            f"(expected ≈ {LOW_BETA_TARGET} for target_beta={LOW_BETA_TARGET})"
        )

    # ------------------------------------------------------------------
    # Future: add test_high_mach_low_beta for the combined extreme once
    # more experience is accumulated with the individual stress tests.
    # ------------------------------------------------------------------
