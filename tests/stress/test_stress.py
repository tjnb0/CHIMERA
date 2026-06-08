"""
Stress tests: verify that the solver behaves correctly under challenging
initial conditions.

Currently tests the GRF/Monte Carlo initialisation by checking that the
realised diagnostics (M_s, beta, gamma) are physically plausible. These
tests will be extended once target_M_s and target_beta are exposed as
CLI arguments to test high-Mach and low-beta regimes explicitly.

The broad bounds used here catch initialisation bugs (e.g. NaN/Inf
parameters, non-positive values) rather than statistical deviations from
the target, which are expected by design since the GRF is random.
"""

import numpy as np
import h5py

GRF_N    = 64
GRF_SEED = 42

M_S_BOUNDS  = (0.01, 20.0)   # sonic Mach number -- physical plausibility
BETA_BOUNDS = (0.01, 100.0)  # plasma beta = 2p/B^2


def _grf_params(run_sim):
    """Run the GRF problem and return the stored parameter dict."""
    h5 = run_sim(problem_type=5, N=GRF_N, seed=GRF_SEED,
                 h5_name=f"stress_grf_N{GRF_N}_s{GRF_SEED}.h5")
    with h5py.File(h5, "r") as f:
        return {
            "gamma": float(f["parameters/gamma"][()]),
            "M_s":   float(f["parameters/M_s"][()]),
            "beta":  float(f["parameters/beta"][()]),
        }


class TestGRFCompletion:
    """GRF runs should complete without crashing at any valid seed."""

    def test_grf_completes_default_seed(self, run_sim):
        h5 = run_sim(problem_type=5, N=GRF_N, seed=GRF_SEED)
        with h5py.File(h5, "r") as f:
            n = len(f["rho"].keys())
        assert n >= 5, f"GRF: only {n} snapshots written"

    def test_grf_completes_different_seed(self, run_sim):
        """Check a second seed to ensure robustness to random variation."""
        h5 = run_sim(problem_type=5, N=GRF_N, seed=GRF_SEED + 1)
        with h5py.File(h5, "r") as f:
            n = len(f["rho"].keys())
        assert n >= 5, f"GRF seed+1: only {n} snapshots written"


class TestGRFParameters:
    """Realised diagnostics must be positive, finite, and physically plausible."""

    def test_gamma_greater_than_one(self, run_sim):
        params = _grf_params(run_sim)
        gamma = params["gamma"]
        assert np.isfinite(gamma), f"gamma is not finite: {gamma}"
        assert gamma > 1.0, f"gamma <= 1: {gamma}"

    def test_mach_number_positive_finite(self, run_sim):
        params = _grf_params(run_sim)
        M_s = params["M_s"]
        assert np.isfinite(M_s), f"M_s is not finite: {M_s}"
        assert M_s > 0.0, f"M_s <= 0: {M_s}"

    def test_mach_number_in_plausible_range(self, run_sim):
        M_s    = _grf_params(run_sim)["M_s"]
        lo, hi = M_S_BOUNDS
        assert lo < M_s < hi, f"M_s = {M_s:.3f} outside [{lo}, {hi}]"

    def test_beta_positive_finite(self, run_sim):
        params = _grf_params(run_sim)
        beta = params["beta"]
        assert np.isfinite(beta), f"beta is not finite: {beta}"
        assert beta > 0.0, f"beta <= 0: {beta}"

    def test_beta_in_plausible_range(self, run_sim):
        beta   = _grf_params(run_sim)["beta"]
        lo, hi = BETA_BOUNDS
        assert lo < beta < hi, f"beta = {beta:.3f} outside [{lo}, {hi}]"

    def test_fields_positive_throughout(self, run_sim):
        """Density and thermal pressure must be positive in all GRF snapshots."""
        h5 = run_sim(problem_type=5, N=GRF_N, seed=GRF_SEED)
        with h5py.File(h5, "r") as f:
            rho_grp = f["rho"]
            p_grp   = f["P"]
            keys    = sorted(rho_grp.keys(),
                             key=lambda k: int(k.replace("SNAPSHOT", "")))
            min_rho = min(np.array(rho_grp[k]).min() for k in keys)
            min_p   = min(np.array(p_grp[k]).min()   for k in keys)
        assert min_rho > 0.0, f"Non-positive density in GRF: min = {min_rho:.3e}"
        assert min_p   > 0.0, f"Non-positive pressure in GRF: min = {min_p:.3e}"

    # ------------------------------------------------------------------
    # Future: add test_high_mach and test_low_beta here once
    # target_M_s and target_beta are exposed as CLI arguments.
    # ------------------------------------------------------------------
