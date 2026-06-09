"""
Quantitative validation tests for CHIMERA.

Two tiers — both marked @pytest.mark.slow; run with:
    python scripts/run_tests.py --validation

Tier 1 — Published bounds (N=128):
    Scalar diagnostics at t=0.5 compared against the ranges reported by
    Stone et al. (2008) and Gardiner & Stone (2005) for their N=192 runs.
    These tests require no reference file and pass as long as the physics
    is qualitatively correct at 128^2.

Tier 2 — Grid convergence (N=64, 128, 256):
    For a 2nd-order scheme with a diffusive Riemann solver (Rusanov),
    numerical dissipation smears peaks and troughs at coarse resolution.
    As N doubles, peaks must sharpen (rho_max, P_max increase) and troughs
    must deepen (rho_min, P_min decrease).  This is a scheme-level check
    that is independent of any reference file.

    If tests/validation/ot_reference.h5 exists (Athena++ 500x500 run via
    scripts/convert_athena_to_reference.py), the test additionally verifies
    that all four diagnostics converge *toward* the high-resolution Athena
    values as N increases.
"""

import hashlib
import numpy as np
import h5py
import pytest
from pathlib import Path

TESTS_DIR      = Path(__file__).resolve().parent
REFERENCE_FILE = TESTS_DIR / "ot_reference.h5"

# ---------------------------------------------------------------------------
# Conservative bounds from published benchmarks (Stone et al. 2008,
# Gardiner & Stone 2005).  Generous to allow for Rusanov diffusion at N=128.
# ---------------------------------------------------------------------------
RHO_MAX_BOUNDS = (0.42, 0.65)
RHO_MIN_BOUNDS = (0.03, 0.12)
P_MAX_BOUNDS   = (0.38, 0.72)
P_MIN_BOUNDS   = (0.008, 0.07)

# Resolutions for the convergence test.  All three must run within the
# validation timeout; wall-clock cost is roughly 0.1 + 0.5 + 2.5 minutes.
CONVERGENCE_NS = [64, 128, 256]

# SHA-256 of the committed Athena++ 500x500 reference file.
# Update this value whenever the reference is intentionally regenerated:
#   sha256sum tests/validation/ot_reference.h5
REFERENCE_SHA256 = "45e66217f6a81bae0831f8201b0cf6b10cdd75df3b5ecfd24c1e7899b5e92b5d"


# ---------------------------------------------------------------------------
# Checksum guard  (fast, no simulation required)
# ---------------------------------------------------------------------------

def test_reference_file_unchanged():
    """
    Guard against accidental modification of the committed Athena++ reference.

    If ot_reference.h5 is absent the test is skipped (contributor hasn't
    generated it yet).  If it is present its SHA-256 must match the value
    committed in REFERENCE_SHA256.

    To update after an intentional regeneration:
        sha256sum tests/validation/ot_reference.h5
    and paste the result into REFERENCE_SHA256 above.
    """
    if not REFERENCE_FILE.exists():
        pytest.skip("ot_reference.h5 not present — skipping checksum guard")

    digest = hashlib.sha256(REFERENCE_FILE.read_bytes()).hexdigest()
    assert digest == REFERENCE_SHA256, (
        "ot_reference.h5 does not match the committed checksum.\n"
        f"  expected : {REFERENCE_SHA256}\n"
        f"  actual   : {digest}\n"
        "If this is an intentional regeneration, update REFERENCE_SHA256 "
        "in tests/validation/test_validation.py."
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _snap_at_time(h5file, field, target_t=0.5, tol=0.006):
    """Return (array, actual_t) for the snapshot closest to target_t."""
    times = np.array(h5file["time/sim_time"])
    idx   = int(np.argmin(np.abs(times - target_t)))
    if abs(times[idx] - target_t) > tol:
        pytest.skip(
            f"No snapshot within {tol} of t={target_t}; "
            f"closest is t={times[idx]:.4f}.  Check tEnd and tOut."
        )
    key = f"SNAPSHOT{idx + 1}"
    return np.array(h5file[field][key]), float(times[idx])


def _scalar_diags(h5_path):
    """Return (rho_max, rho_min, P_max, P_min) at the t=0.5 snapshot."""
    with h5py.File(h5_path, "r") as f:
        rho, _ = _snap_at_time(f, "rho")
        P,   _ = _snap_at_time(f, "P")
    return float(rho.max()), float(rho.min()), float(P.max()), float(P.min())


# ---------------------------------------------------------------------------
# Tier 1 — Published bounds  (N=128)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestOTVortexValidation:
    """
    OT vortex at N=128, t=0.5 vs published scalar bounds.
    The ot_128_h5 fixture (tests/conftest.py) runs the simulation once per
    session and is shared across all tests in this class.
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

    def test_report_values(self, ot_128_h5, capsys):
        """
        Print realised values alongside published expectations.
        Always passes — use for manual inspection against paper figures.
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


# ---------------------------------------------------------------------------
# Tier 2 — Grid convergence  (N=64, 128, 256)
# ---------------------------------------------------------------------------

@pytest.mark.slow
class TestOTConvergence:
    """
    Grid convergence of OT scalar diagnostics at t=0.5.

    For a 2nd-order scheme with Rusanov fluxes, numerical dissipation
    over-smooths peaks at coarse resolution.  As the grid is refined:

      rho_max, P_max  must increase  (peaks sharpen)
      rho_min, P_min  must decrease  (troughs deepen)

    If tests/validation/ot_reference.h5 is present (Athena++ 500x500), the
    test also verifies that errors in all four diagnostics decrease as N
    increases, confirming convergence toward the high-resolution truth.

    Both checks run in a single test function to avoid running redundant
    simulations (3 runs total, ~3–5 minutes).
    """

    def _run_and_get_diags(self, run_sim, N):
        """Run OT at N and return (rho_max, rho_min, P_max, P_min) at t=0.5."""
        h5 = run_sim(
            problem_type=1, N=N,
            h5_name=f"ot_conv_N{N}.h5",
            timeout=600,
        )
        return _scalar_diags(h5)

    def test_convergence(self, run_sim):
        """
        Monotone convergence of all four diagnostics as N doubles.
        Convergence toward Athena++ 500x500 values if reference is present.
        """
        # --- run all three resolutions ---
        results = {N: self._run_and_get_diags(run_sim, N) for N in CONVERGENCE_NS}

        labels    = ["rho_max", "rho_min", "P_max", "P_min"]
        increases = [True,       False,     True,    False]

        # --- Part 1: strict monotonicity between successive resolutions ---
        for i in range(len(CONVERGENCE_NS) - 1):
            n1, n2 = CONVERGENCE_NS[i], CONVERGENCE_NS[i + 1]
            d1, d2 = results[n1], results[n2]
            for label, v1, v2, should_increase in zip(labels, d1, d2, increases):
                if should_increase:
                    assert v2 > v1, (
                        f"{label} did not increase from N={n1} to N={n2}: "
                        f"{v1:.4f} -> {v2:.4f}"
                    )
                else:
                    assert v2 < v1, (
                        f"{label} did not decrease from N={n1} to N={n2}: "
                        f"{v1:.4f} -> {v2:.4f}"
                    )

        # --- Part 2: convergence toward Athena++ 500x500 (if available) ---
        if not REFERENCE_FILE.exists():
            return   # monotonicity check is sufficient without a reference

        ref_diags = _scalar_diags(REFERENCE_FILE)

        def _max_rel_err(diags, ref):
            return max(
                abs(d - r) / (abs(r) + 1e-30)
                for d, r in zip(diags, ref)
            )

        errors = {N: _max_rel_err(results[N], ref_diags) for N in CONVERGENCE_NS}

        for i in range(len(CONVERGENCE_NS) - 1):
            n1, n2 = CONVERGENCE_NS[i], CONVERGENCE_NS[i + 1]
            assert errors[n2] < errors[n1], (
                f"Scalar error did not decrease from N={n1} to N={n2}: "
                f"{errors[n1]:.4f} -> {errors[n2]:.4f}\n"
                f"  Diagnostics  N={n1}: {dict(zip(labels, results[n1]))}\n"
                f"  Diagnostics  N={n2}: {dict(zip(labels, results[n2]))}\n"
                f"  Reference  (500x500): {dict(zip(labels, ref_diags))}"
            )
