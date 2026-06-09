"""
Structural tests: verify that HDF5 output has the correct format,
field names, dimensions, and time ordering.

These tests check that the solver produces well-formed output regardless
of whether the physics is correct. They run at low resolution (N=32)
and complete in seconds.
"""

import numpy as np
import h5py

FIELDS  = ["rho", "P", "Bx", "By", "Vx", "Vy"]
N_TEST  = 32
PARAMS  = ["parameters/gamma", "parameters/M_s", "parameters/beta"]


class TestHDF5Structure:
    """Check that every required dataset is present and well-formed."""

    def test_all_fields_present(self, run_sim):
        h5 = run_sim(problem_type=1, N=N_TEST)
        with h5py.File(h5, "r") as f:
            for field in FIELDS:
                assert field in f, f"Missing field group: {field}"

    def test_time_dataset_present(self, run_sim):
        h5 = run_sim(problem_type=1, N=N_TEST)
        with h5py.File(h5, "r") as f:
            assert "time/sim_time" in f, "Missing time/sim_time"

    def test_physics_parameters_present(self, run_sim):
        """gamma, M_s, and beta must be written for all problem types."""
        h5 = run_sim(problem_type=5, N=N_TEST, seed=1)
        with h5py.File(h5, "r") as f:
            for param in PARAMS:
                assert param in f, f"Missing parameter: {param}"

    def test_snapshot_shape_matches_N(self, run_sim):
        """Every snapshot must be exactly (N, N)."""
        h5 = run_sim(problem_type=1, N=N_TEST)
        with h5py.File(h5, "r") as f:
            grp  = f["rho"]
            key  = sorted(grp.keys(),
                          key=lambda k: int(k.replace("SNAPSHOT", "")))[0]
            shape = np.array(grp[key]).shape
        assert shape == (N_TEST, N_TEST), \
            f"Expected ({N_TEST},{N_TEST}), got {shape}"

    def test_snapshot_counts_consistent_across_fields(self, run_sim):
        """All six fields must have the same number of snapshots."""
        h5 = run_sim(problem_type=1, N=N_TEST)
        with h5py.File(h5, "r") as f:
            counts = {field: len(f[field].keys()) for field in FIELDS}
        values = list(counts.values())
        assert len(set(values)) == 1, \
            f"Inconsistent snapshot counts: {counts}"

    def test_time_array_length_matches_snapshots(self, run_sim):
        """Length of time/sim_time must equal number of rho snapshots."""
        h5 = run_sim(problem_type=1, N=N_TEST)
        with h5py.File(h5, "r") as f:
            n_snaps = len(f["rho"].keys())
            n_times = len(np.array(f["time/sim_time"]))
        assert n_snaps == n_times, \
            f"Snapshot count ({n_snaps}) != time array length ({n_times})"

    def test_times_monotonically_increasing(self, run_sim):
        """Snapshot times must be strictly increasing."""
        h5 = run_sim(problem_type=1, N=N_TEST)
        with h5py.File(h5, "r") as f:
            t = np.array(f["time/sim_time"])
        assert np.all(np.diff(t) > 0), \
            f"Times not monotonically increasing: {t[:5]} ..."

    def test_at_least_five_snapshots_written(self, run_sim):
        """Confirms the simulation ran meaningfully, not just one step."""
        h5 = run_sim(problem_type=1, N=N_TEST)
        with h5py.File(h5, "r") as f:
            n = len(f["rho"].keys())
        assert n >= 5, f"Only {n} snapshots written"

    def test_p_dataset_is_thermal_pressure(self, run_sim):
        """
        The P dataset must store thermal pressure p = P_total - 0.5*B^2,
        not total pressure P*.

        Discrimination method: for the OT vortex the initial thermal
        pressure is uniform at p_0 = 5/(12π) ≈ 0.133.  The mean
        magnetic pressure over the periodic domain is
        mean(B^2/2) = 1/(8π) ≈ 0.040, so mean total P* ≈ 0.173.

        At the first snapshot (t ≈ tOut = 0.01), the flow has barely
        evolved, so if P stores thermal pressure its domain mean should
        still be close to 0.133 (<< 0.16).  If P accidentally stored
        total pressure, the mean would be ~0.173 (>> 0.16).
        """
        h5 = run_sim(problem_type=1, N=N_TEST)
        with h5py.File(h5, "r") as f:
            keys    = sorted(f["P"].keys(),
                             key=lambda k: int(k.replace("SNAPSHOT", "")))
            P_first = np.array(f["P"][keys[0]])
            P_last  = np.array(f["P"][keys[-1]])

        # Positivity (applies at all times)
        assert P_first.min() > 0.0, \
            f"Thermal pressure not positive at first snapshot: min = {P_first.min():.3e}"
        assert P_last.min() > 0.0, \
            f"Thermal pressure not positive at last snapshot: min = {P_last.min():.3e}"

        # Domain-mean anchor: separates thermal from total pressure.
        # Threshold 0.16 lies between p_0 ≈ 0.133 (thermal) and P*_0 ≈ 0.173 (total).
        p_thermal_init = 5.0 / (12.0 * np.pi)   # ≈ 0.133
        THRESHOLD      = 0.16
        assert P_first.mean() < THRESHOLD, (
            f"First-snapshot mean P = {P_first.mean():.4f} >= {THRESHOLD}. "
            f"Expected near {p_thermal_init:.4f} for thermal pressure; "
            f"mean total pressure would be ~0.173 — P may include B²/2."
        )
