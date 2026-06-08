"""
Quantitative validation tests for CHIMERA.

These tests confirm that the solver produces results consistent with
published benchmarks, not just that it runs without crashing.

Marked @pytest.mark.slow because they require the N=128 OT vortex run
(~10-30 seconds). Skipped by default; run with:
    pytest tests/validation/ -m slow
    python scripts/run_tests.py --validation

Reference data (tests/validation/ot_reference.h5):
    Generated once by  python scripts/generate_ot_reference.py  after
    visually validating the N=128 output against the published figures.
    Commit the file so all future runs compare against the same baseline.

Published scalar bounds (Stone et al. 2008, Gardiner & Stone 2005):
    At t=0.5, N=128, gamma=5/3:
      rho_max ~ 0.52,  rho_min ~ 0.07
      p_max   ~ 0.53,  p_min   ~ 0.02
"""

import numpy as np
import h5py
import pytest
from pathlib import Path

TESTS_DIR      = Path(__file__).resolve().parent
REFERENCE_FILE = TESTS_DIR / "ot_reference.h5"

# Conservative bounds from published benchmarks (generous to allow for
# minor scheme differences and N=128 resolution effects).
RHO_MAX_BOUNDS = (0.42, 0.65)
RHO_MIN_BOUNDS = (0.03, 0.12)
P_MAX_BOUNDS   = (0.38, 0.72)
P_MIN_BOUNDS   = (0.008, 0.07)

# Tolerances for comparison against the stored reference
L2_DENSITY_TOL = 0.03   # 3% relative L2 error
L1_PROFILE_TOL = 0.05   # 5% mean relative L1 error on pressure slice


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _snap_at_time(h5file, field, target_t=0.5, tol=0.006):
    """Return the snapshot array of `field` closest to target_t."""
    times = np.array(h5file["time/sim_time"])
    idx   = int(np.argmin(np.abs(times - target_t)))
    if abs(times[idx] - target_t) > tol:
        pytest.skip(
            f"No snapshot within {tol} of t={target_t}; "
            f"closest is t={times[idx]:.4f}. Check tEnd and tOut."
        )
    key = f"SNAPSHOT{idx + 1}"
    return np.array(h5file[field][key]), float(times[idx])


# ---------------------------------------------------------------------------
# Validation tests
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestOTVortexValidation:
    """
    OT vortex at N=128, t=0.5 vs published bounds and stored reference.
    The ot_128_h5 fixture (defined in tests/conftest.py) runs the simulation
    exactly once per session and is shared across all tests in this class.
    """

    def test_density_max_in_published_bounds(self, ot_128_h5):
        """max(rho) must fall within bounds from Stone et al. (2008)."""
        with h5py.File(ot_128_h5, "r") as f:
            rho, t = _snap_at_time(f, "rho")
        lo, hi = RHO_MAX_BOUNDS
        assert lo < rho.max() < hi, (
            f"max(rho) = {rho.max():.4f} outside [{lo}, {hi}] at t={t:.4f}"
        )

    def test_density_min_in_published_bounds(self, ot_128_h5):
        """min(rho) must fall within bounds from Stone et al. (2008)."""
        with h5py.File(ot_128_h5, "r") as f:
            rho, t = _snap_at_time(f, "rho")
        lo, hi = RHO_MIN_BOUNDS
        assert lo < rho.min() < hi, (
            f"min(rho) = {rho.min():.4f} outside [{lo}, {hi}] at t={t:.4f}"
        )

    def test_pressure_max_in_published_bounds(self, ot_128_h5):
        """max(p_thermal) must fall within bounds from Stone et al. (2008)."""
        with h5py.File(ot_128_h5, "r") as f:
            P, t = _snap_at_time(f, "P")
        lo, hi = P_MAX_BOUNDS
        assert lo < P.max() < hi, (
            f"max(P) = {P.max():.4f} outside [{lo}, {hi}] at t={t:.4f}"
        )

    def test_pressure_min_in_published_bounds(self, ot_128_h5):
        """min(p_thermal) must fall within bounds from Stone et al. (2008)."""
        with h5py.File(ot_128_h5, "r") as f:
            P, t = _snap_at_time(f, "P")
        lo, hi = P_MIN_BOUNDS
        assert lo < P.min() < hi, (
            f"min(P) = {P.min():.4f} outside [{lo}, {hi}] at t={t:.4f}"
        )

    def test_density_l2_error_vs_reference(self, ot_128_h5):
        """
        Relative L2 error of density vs stored reference must be within
        L2_DENSITY_TOL. Skipped if ot_reference.h5 does not exist.

        Generate the reference with:
            python scripts/generate_ot_reference.py
        """
        if not REFERENCE_FILE.exists():
            pytest.skip(
                f"Reference file not found: {REFERENCE_FILE}\n"
                f"Run: python scripts/generate_ot_reference.py"
            )
        with h5py.File(REFERENCE_FILE, "r") as ref:
            rho_ref, _ = _snap_at_time(ref, "rho")
        with h5py.File(ot_128_h5, "r") as f:
            rho, t = _snap_at_time(f, "rho")
        rms_ref = np.sqrt(np.mean(rho_ref**2))
        l2_err  = np.sqrt(np.mean((rho - rho_ref)**2)) / (rms_ref + 1e-30)
        assert l2_err < L2_DENSITY_TOL, (
            f"Density L2 error vs reference = {l2_err:.4f} "
            f"(threshold {L2_DENSITY_TOL}) at t={t:.4f}"
        )

    def test_pressure_profile_vs_reference(self, ot_128_h5):
        """
        Horizontal pressure profile at y=0.5 (middle row) vs stored reference.
        Skipped if ot_reference.h5 does not exist.
        """
        if not REFERENCE_FILE.exists():
            pytest.skip(
                f"Reference file not found: {REFERENCE_FILE}\n"
                f"Run: python scripts/generate_ot_reference.py"
            )
        with h5py.File(REFERENCE_FILE, "r") as ref:
            P_ref, _ = _snap_at_time(ref, "P")
        with h5py.File(ot_128_h5, "r") as f:
            P, t = _snap_at_time(f, "P")

        N = P.shape[0]
        profile     = P[:, N // 2]
        profile_ref = P_ref[:, N // 2]
        mean_ref    = np.mean(np.abs(profile_ref)) + 1e-30
        l1_err      = np.mean(np.abs(profile - profile_ref)) / mean_ref
        assert l1_err < L1_PROFILE_TOL, (
            f"Pressure profile L1 error vs reference = {l1_err:.4f} "
            f"(threshold {L1_PROFILE_TOL}) at t={t:.4f}"
        )

    def test_report_values(self, ot_128_h5, capsys):
        """
        Print actual values alongside published expectations. Always passes.
        Use this for manual inspection against paper figures.
        """
        with h5py.File(ot_128_h5, "r") as f:
            rho, t = _snap_at_time(f, "rho")
            P,   _ = _snap_at_time(f, "P")
        print(
            f"\n  OT vortex at t={t:.4f}, N=128:\n"
            f"    rho: min={rho.min():.4f}  max={rho.max():.4f}\n"
            f"    P:   min={P.min():.4f}  max={P.max():.4f}\n"
            f"\n  Expected (Stone et al. 2008, Gardiner & Stone 2005):\n"
            f"    rho: min~0.07  max~0.52\n"
            f"    P:   min~0.02  max~0.53"
        )
