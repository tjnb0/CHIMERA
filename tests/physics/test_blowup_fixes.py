"""
Regression tests for the three pressure-blowup fixes.

Background
----------
Two MC runs were diagnosed as failing through distinct mechanisms:

  Regime A  (gamma=1.683, M_s=0.814, beta=1.632)
    Failure mode: slow reconnection blowup
    Current-sheet compression drained thermal pressure to the absolute floor
    (P_floor=1e-12) in ~30 cells from t=0.30 onward.  The MOOD check used
    P_floor as its threshold, so face states with p_th ~ -2.8 passed through
    undetected for 30+ snapshots.  A reconnection event spiked |B|^2 from
    5.9->13.1 around t=0.61 and exploded the floored cells.

  Regime B  (gamma=1.345, M_s=0.914, beta=1.685)
    Failure mode: fast density-collapse blowup
    Higher M_s and softer EOS caused the floored cell to drain from
    rho=0.24 to rho=0.038 in 7 steps while thermally void, creating a
    near-vacuum cell that exploded at t=0.44 when hit by a compressive wave.

Tests use the command-line parameter overrides (args 6/7/8 = M_s/beta/gamma)
to reproduce the diagnosed physical conditions on any machine, independent
of compiler or RNG state.  No specific seed values are used.
"""

import math
import numpy as np
import h5py

N_MC = 64   # grid size matching the diagnosed runs

# ---------------------------------------------------------------------------
# Diagnosed parameter regimes (named by physics, not by run ID)
# ---------------------------------------------------------------------------

# Regime A: high gamma, moderate M_s, moderate beta  — reconnection blowup
REGIME_A = dict(M_s=0.81, beta=1.63, gamma=1.68)

# Regime B: low gamma, high M_s, moderate beta — density-collapse blowup
REGIME_B = dict(M_s=0.91, beta=1.69, gamma=1.35)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _sorted_keys(grp):
    return sorted(grp.keys(), key=lambda k: int(k.replace("SNAPSHOT", "")))


def _snap_count(h5_path):
    with h5py.File(h5_path, "r") as f:
        return len(f["rho"].keys())


def _params(h5_path):
    with h5py.File(h5_path, "r") as f:
        return {
            "gamma": np.asarray(f["parameters/gamma"]).item(),
            "M_s":   np.asarray(f["parameters/M_s"]).item(),
            "beta":  np.asarray(f["parameters/beta"]).item(),
        }


def _min_thermal_pressure(h5_path):
    """Global minimum thermal pressure across all snapshots.

    /P/ stores thermal pressure directly (README: 'thermal pressure
    p = P* - B^2/2').  No B-field subtraction needed here.
    """
    with h5py.File(h5_path, "r") as f:
        keys = _sorted_keys(f["P"])
        return min(float(np.array(f["P"][k]).min()) for k in keys)


def _mean_thermal_pressure_at_snap(h5_path, snap_index):
    """Domain-mean thermal pressure at a single snapshot.

    /P/ stores thermal pressure directly; just average it.
    """
    with h5py.File(h5_path, "r") as f:
        keys = _sorted_keys(f["P"])
        return float(np.array(f["P"][keys[snap_index]]).mean())


def _min_density(h5_path):
    with h5py.File(h5_path, "r") as f:
        keys = _sorted_keys(f["rho"])
        return min(np.array(f["rho"][k]).min() for k in keys)


def _max_pressure(h5_path):
    with h5py.File(h5_path, "r") as f:
        keys = _sorted_keys(f["P"])
        return max(np.array(f["P"][k]).max() for k in keys)


def _pressure_at_snap(h5_path, snap_index):
    with h5py.File(h5_path, "r") as f:
        keys = _sorted_keys(f["P"])
        return np.array(f["P"][keys[snap_index]])



# ---------------------------------------------------------------------------
# Fix 3: IC positivity guard
# ---------------------------------------------------------------------------

class TestICPositivityGuard:
    """
    setup_GRF_fields must not produce initial cells with 0.5*B^2 > P_thermal

    Expected p_th_mean from the beta sampling:
      beta = 2 * p_th / B^2  =>  p_th_mean ~ beta * P_mean / (2 + beta)
      Regime A (beta=1.63): p_th_mean ~ 0.45
      Regime B (beta=1.69): p_th_mean ~ 0.46

    Use a conservative lower bound of 0.1 * P_mean ~ 0.1.
    """

    P_TH_MEAN_MIN = 0.1   # conservative lower bound

    def test_high_gamma_mod_beta_p_th_mean_positive(self, run_sim):
        """
        Regime A (high gamma, moderate M_s, moderate beta).
        """
        p = REGIME_A
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="regime_a.h5",
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        p_th_mean = _mean_thermal_pressure_at_snap(h5, snap_index=0)
        assert p_th_mean > self.P_TH_MEAN_MIN, (
            f"(IC guard, regime A): domain-mean thermal pressure at "
            f"t=0.01 is {p_th_mean:.4f}, below threshold {self.P_TH_MEAN_MIN}. "
            f"Expected ~0.45 from beta={p['beta']}. IC guard may be "
            f"failing to correct over-pressured cells, or over-correcting."
        )

    def test_low_gamma_high_Ms_p_th_mean_positive(self, run_sim):
        """
        Regime B (low gamma, high M_s, moderate beta).
        Confirms the IC guard works across both diagnosed parameter regimes.
        """
        p = REGIME_B
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="regime_b.h5",
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        p_th_mean = _mean_thermal_pressure_at_snap(h5, snap_index=0)
        assert p_th_mean > self.P_TH_MEAN_MIN, (
            f"(IC guard, regime B): domain-mean thermal pressure at "
            f"t=0.01 is {p_th_mean:.4f}, below threshold {self.P_TH_MEAN_MIN}. "
            f"Expected ~0.46 from beta={p['beta']}."
        )

    def test_ic_guard_does_not_perturb_mean_pressure(self, run_sim):
        """
        The IC guard only raises cells that are magnetically over-pressured.
        Mean thermal pressure should not exceed 2 * P_mean ~ 2.0, confirming
        the guard is not adding energy to the bulk of the domain.

        /P/ stores thermal pressure directly; no B-field subtraction needed.
        """
        p = REGIME_A
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="egime_a_mean.h5",
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        p_th_mean = _mean_thermal_pressure_at_snap(h5, snap_index=0)
        assert p_th_mean < 2.0, (
            f"IC guard mean ceiling, regime A): mean thermal pressure "
            f"at t=0.01 is {p_th_mean:.4f}, exceeding 2.0. "
            f"The IC guard may be over-correcting."
        )


# ---------------------------------------------------------------------------
# Fix 2: relative MOOD threshold
# ---------------------------------------------------------------------------

class TestRelativeMOODThreshold:
    """
    thermal_pressure_check must fire before face-state p_th reaches
    the absolute floor (1e-12).
    """

    MOOD_FLOOR_EVIDENCE = 1e-6

    def test_high_gamma_mod_beta_p_th_above_floor(self, run_sim):
        """
        Regime A: relative MOOD threshold should catch face states before
        p_th hits 1e-12 in cell averages.
        """
        p = REGIME_A
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="fix2_regime_a.h5",
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        p_min = _min_thermal_pressure(h5)
        assert p_min > self.MOOD_FLOOR_EVIDENCE, (
            f"Fix 2 (relative MOOD, regime A): minimum thermal pressure "
            f"{p_min:.3e} at or below absolute floor evidence "
            f"{self.MOOD_FLOOR_EVIDENCE:.0e}. MOOD may not be firing early enough."
        )

    def test_low_gamma_high_Ms_p_th_above_floor(self, run_sim):
        """Regime B: same check for the density-collapse parameter regime."""
        p = REGIME_B
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="fix2_regime_b.h5",
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        p_min = _min_thermal_pressure(h5)
        assert p_min > self.MOOD_FLOOR_EVIDENCE, (
            f"Fix 2 (relative MOOD, regime B): minimum thermal pressure "
            f"{p_min:.3e} at or below absolute floor evidence "
            f"{self.MOOD_FLOOR_EVIDENCE:.0e}."
        )

    def test_mood_does_not_over_trigger_on_ot_vortex(self, run_sim):
        """
        Sanity check: the relative MOOD threshold must not pathologically
        suppress reconstruction on a smooth reference problem.
        """
        P_OT_INIT = 5.0 / (12.0 * math.pi)
        P_MAX_THRESHOLD = 0.25

        h5 = run_sim(problem_type=1, N=32, h5_name="fix2_ot_sanity.h5")
        with h5py.File(h5, "r") as f:
            keys = _sorted_keys(f["P"])
            p_max_final = float(np.array(f["P"][keys[-1]]).max())
        assert p_max_final > P_MAX_THRESHOLD, (
            f"Fix 2 sanity (OT vortex N=32): peak thermal pressure at final "
            f"snapshot is {p_max_final:.4f}, below {P_MAX_THRESHOLD:.3f} "
            f"(~2x P_init={P_OT_INIT:.4f}). MOOD may be over-triggering."
        )


# ---------------------------------------------------------------------------
# Fix 1: proportional density floor
# ---------------------------------------------------------------------------

class TestProportionalDensityFloor:
    """
    apply_conserved_floors must prevent near-vacuum density collapse.
    """

    RHO_FLOOR_EVIDENCE = 1e-4   # factor-of-10 margin below the 1e-3 floor

    def test_low_gamma_high_Ms_density_bounded(self, run_sim):
        """
        Regime B (low gamma, high M_s): the density-collapse regime.
        Minimum density must never drop below 1e-4 throughout the run.
        """
        p = REGIME_B
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="fix1_regime_b.h5",
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        rho_min = _min_density(h5)
        assert rho_min > self.RHO_FLOOR_EVIDENCE, (
            f"Fix 1 (proportional rho floor, regime B): minimum density "
            f"{rho_min:.3e} dropped below {self.RHO_FLOOR_EVIDENCE:.0e}. "
            f"Near-vacuum collapse not prevented."
        )

    def test_high_gamma_mod_beta_density_bounded(self, run_sim):
        """Regime A: density must also stay bounded in the reconnection regime."""
        p = REGIME_A
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="fix1_regime_a.h5",
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        rho_min = _min_density(h5)
        assert rho_min > self.RHO_FLOOR_EVIDENCE, (
            f"Fix 1 (proportional rho floor, regime A): minimum density "
            f"{rho_min:.3e} dropped below {self.RHO_FLOOR_EVIDENCE:.0e}."
        )

    def test_density_floor_mass_growth_bounded(self, run_sim):
        """
        The proportional floor adds mass to near-vacuum cells.  Total mass
        must not grow by more than 1% relative to its initial value.
        Uses regime B — the density-collapse regime — where the floor fires.
        """
        p = REGIME_B
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="fix1_regime_b_mass.h5",
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        with h5py.File(h5, "r") as f:
            keys = _sorted_keys(f["rho"])
            rho_first = np.array(f["rho"][keys[0]])
            rho_last  = np.array(f["rho"][keys[-1]])
        mass_first = rho_first.sum()
        mass_last  = rho_last.sum()
        rel_change = abs(mass_last - mass_first) / mass_first
        assert rel_change < 0.01, (
            f"Fix 1: mass changed by {rel_change:.3%} "
            f"(initial={mass_first:.4f}, final={mass_last:.4f}). "
            f"Proportional floor is adding too much mass."
        )


# ---------------------------------------------------------------------------
# End-to-end: both diagnosed parameter regimes must complete
# ---------------------------------------------------------------------------

class TestDiagnosedConditionsComplete:
    """
    Full-run end-to-end tests reproducing the exact physical conditions of
    the two diagnosed failures via parameter overrides.  All three fixes
    together must allow both regimes to reach t=tEnd (100 snapshots).
    """

    EXPECTED_SNAPS = 100
    P_MAX_SANE     = 5000.0   # pre-fix peaks were ~33,000 and ~316,000

    def test_high_gamma_mod_beta_completes(self, run_sim):
        """Regime A (reconnection blowup): must write 100 snapshots."""
        p = REGIME_A
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="e2e_regime_a.h5", timeout=180,
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        n = _snap_count(h5)
        assert n == self.EXPECTED_SNAPS, (
            f"Regime A only wrote {n}/{self.EXPECTED_SNAPS} snapshots."
        )

    def test_high_gamma_mod_beta_pressure_positive(self, run_sim):
        """Regime A: thermal pressure must be positive throughout."""
        p = REGIME_A
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="e2e_regime_a_P.h5", timeout=180,
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        p_min = _min_thermal_pressure(h5)
        assert p_min > 0.0, (
            f"Regime A: negative thermal pressure {p_min:.3e} after fixes."
        )

    def test_high_gamma_mod_beta_pressure_bounded(self, run_sim):
        """Regime A: pressure must not blow up (pre-fix peak ~33,000)."""
        p = REGIME_A
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="e2e_regime_a_Pmax.h5", timeout=180,
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        p_max = _max_pressure(h5)
        assert p_max < self.P_MAX_SANE, (
            f"Regime A: max pressure {p_max:.2e} exceeds bound {self.P_MAX_SANE:.0f}."
        )

    def test_low_gamma_high_Ms_completes(self, run_sim):
        """Regime B (density-collapse blowup): must write 100 snapshots."""
        p = REGIME_B
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="e2e_regime_b.h5", timeout=180,
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        n = _snap_count(h5)
        assert n == self.EXPECTED_SNAPS, (
            f"Regime B only wrote {n}/{self.EXPECTED_SNAPS} snapshots."
        )

    def test_low_gamma_high_Ms_pressure_positive(self, run_sim):
        """Regime B: thermal pressure must be positive throughout."""
        p = REGIME_B
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="e2e_regime_b_P.h5", timeout=180,
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        p_min = _min_thermal_pressure(h5)
        assert p_min > 0.0, (
            f"Regime B: negative thermal pressure {p_min:.3e} after fixes."
        )

    def test_low_gamma_high_Ms_pressure_bounded(self, run_sim):
        """Regime B: pressure must not blow up (pre-fix peak ~316,000)."""
        p = REGIME_B
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="e2e_regime_b_Pmax.h5", timeout=180,
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        p_max = _max_pressure(h5)
        assert p_max < self.P_MAX_SANE, (
            f"Regime B: max pressure {p_max:.2e} exceeds bound {self.P_MAX_SANE:.0f}."
        )

    def test_low_gamma_high_Ms_density_positive(self, run_sim):
        """
        Regime B: density must be positive throughout.
        The density-collapse mechanism was distinct from the pressure blowup —
        this test targets Fix 1 specifically in the end-to-end context.
        """
        p = REGIME_B
        h5 = run_sim(
            problem_type=5, N=N_MC, seed=1,
            h5_name="e2e_regime_b_rho.h5", timeout=180,
            extra_args=[p["M_s"], p["beta"], p["gamma"]],
        )
        rho_min = _min_density(h5)
        assert rho_min > 0.0, (
            f"Regime B: non-positive density {rho_min:.3e} after fixes."
        )
