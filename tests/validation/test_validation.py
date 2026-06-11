"""
Quantitative validation tests for CHIMERA.

Two tiers - both marked @pytest.mark.slow; run with:
    python scripts/run_tests.py --validation

Tier 1 - Published bounds (N=128):
    Scalar diagnostics at t=0.5 compared against the ranges reported by
    Stone et al. (2008) and Gardiner & Stone (2005) for their N=192 runs.
    These tests require no reference file and pass as long as the physics
    is qualitatively correct at 128^2.

Tier 2 - Grid convergence (N=64, 128, 256):
    Self-convergence: use N=256 as the high-resolution reference.  Coarsen
    it to N=64 and N=128 by block-averaging and verify that the L2 error
    decreases as N doubles with convergence rate > 1.  This check is
    scheme-agnostic — it makes no assumption about which Riemann solver is
    active.

    If tests/validation/ot_reference.h5 exists (Athena++ 500x500 run via
    scripts/convert_athena_to_reference.py), the test additionally verifies
    that L2 errors against the external reference decrease as N increases.
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
# Gardiner & Stone 2005).  Generous to allow for numerical diffusion at N=128.
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


def _coarsen(arr, factor):
    """
    Block-average a 2D array by an integer factor in both dimensions.

    Example: (256, 256) with factor=4 -> (64, 64).
    Requires arr.shape[0] % factor == 0.
    """
    n = arr.shape[0]
    m = n // factor
    return arr.reshape(m, factor, m, factor).mean(axis=(1, 3))


def _l2_error(a, b):
    """RMS element-wise difference between two same-shape arrays."""
    return float(np.sqrt(np.mean((a - b) ** 2)))


def _interp_to_grid(arr, n_dst):
    """
    Bilinear interpolation of arr (n_src x n_src) onto an n_dst x n_dst
    uniform cell-centred grid covering [0, 1]^2.

    Used to compare against the 500x500 Athena++ reference, where 500 is
    not an integer multiple of 64 / 128 / 256.
    """
    n_src = arr.shape[0]
    x_src = (np.arange(n_src) + 0.5) / n_src   # source cell centres
    x_dst = (np.arange(n_dst) + 0.5) / n_dst   # destination cell centres

    # Interpolate along x-axis (axis=1) for every row
    tmp = np.empty((n_src, n_dst))
    for i in range(n_src):
        tmp[i] = np.interp(x_dst, x_src, arr[i])

    # Interpolate along y-axis (axis=0) for every column
    result = np.empty((n_dst, n_dst))
    for j in range(n_dst):
        result[:, j] = np.interp(x_dst, x_src, tmp[:, j])

    return result


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
    L2 grid convergence of the OT vortex at t=0.5.

    Part 1 — Self-convergence (no reference file required):
        Uses N=256 as the high-resolution surrogate truth.  The N=256
        solution is coarsened to N=64 and N=128 by block-averaging (factor
        4 and 2 respectively, both exact integer ratios).  The L2 error of
        each coarser run against the coarsened N=256 solution must decrease
        as N doubles, with rate = log2(L2_64 / L2_128) > 1.0.

        Block-averaging is used rather than subsampling so that the coarsened
        reference represents the same cell-averaged quantity as the coarser
        simulation, minimising aliasing error in the rate estimate.

        This check is scheme-agnostic: it passes for any correct solver,
        regardless of how diffusive it is.

    Part 2 — Convergence toward Athena++ 500x500 (if ot_reference.h5 exists):
        The reference is interpolated to each simulation grid via bilinear
        interpolation (500 is not divisible by 64/128/256, so block-averaging
        cannot be used here).  L2 errors against the external reference must
        decrease as N increases.
    """

    def test_convergence(self, run_sim):
        """
        L2 convergence in rho and thermal pressure as N doubles.

        Runs OT at N=64, 128, 256 (total ~3–5 min).
        Asserts self-convergence rate > 1.0 for both fields,
        then (if ot_reference.h5 is present) asserts convergence
        toward the Athena++ 500x500 reference.
        """
        # --- Run all three resolutions ---
        h5s = {
            N: run_sim(problem_type=1, N=N, h5_name=f"ot_conv_N{N}.h5", timeout=600)
            for N in CONVERGENCE_NS
        }

        # Load rho and P arrays at t≈0.5 for each resolution
        fields = {}
        for N, h5 in h5s.items():
            with h5py.File(h5, "r") as f:
                rho, _ = _snap_at_time(f, "rho")
                P,   _ = _snap_at_time(f, "P")
            fields[N] = {"rho": rho, "P": P}

        # --- Part 1: Self-convergence using N=256 as reference ---
        # Block-average N=256 to each coarser resolution and compute L2 error.
        # Convergence rate = log2(L2_64 / L2_128) must exceed 1.0.
        ref_N = CONVERGENCE_NS[-1]   # 256
        for fname in ("rho", "P"):
            ref_arr = fields[ref_N][fname]
            errors  = {
                N: _l2_error(fields[N][fname], _coarsen(ref_arr, ref_N // N))
                for N in CONVERGENCE_NS[:-1]   # 64, 128
            }

            n1, n2 = CONVERGENCE_NS[0], CONVERGENCE_NS[1]   # 64, 128
            rate   = np.log2(errors[n1] / errors[n2])

            assert errors[n2] < errors[n1], (
                f"L2({fname}) did not decrease from N={n1} to N={n2}: "
                f"{errors[n1]:.3e} -> {errors[n2]:.3e}"
            )
            assert rate > 1.0, (
                f"L2({fname}) convergence rate too low: {rate:.2f} (expected > 1.0)\n"
                f"  N={n1}: L2={errors[n1]:.3e},  N={n2}: L2={errors[n2]:.3e}"
            )

        # --- Part 2: Convergence toward Athena++ 500x500 (if available) ---
        if not REFERENCE_FILE.exists():
            return

        ref_fields = {}
        with h5py.File(REFERENCE_FILE, "r") as f:
            ref_fields["rho"], _ = _snap_at_time(f, "rho")
            ref_fields["P"],   _ = _snap_at_time(f, "P")

        for fname in ("rho", "P"):
            errors = {
                N: _l2_error(fields[N][fname], _interp_to_grid(ref_fields[fname], N))
                for N in CONVERGENCE_NS
            }

            for i in range(len(CONVERGENCE_NS) - 1):
                n1, n2 = CONVERGENCE_NS[i], CONVERGENCE_NS[i + 1]
                assert errors[n2] < errors[n1], (
                    f"Athena++ L2({fname}) did not decrease from N={n1} to N={n2}: "
                    f"{errors[n1]:.3e} -> {errors[n2]:.3e}"
                )
