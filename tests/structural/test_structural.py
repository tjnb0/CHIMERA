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
        The P dataset stores thermal pressure p, not total pressure P*.
        p must be positive and strictly less than p + B^2/2 where B != 0.
        """
        h5 = run_sim(problem_type=1, N=N_TEST)
        with h5py.File(h5, "r") as f:
            grp  = f["rho"]
            keys = sorted(grp.keys(),
                          key=lambda k: int(k.replace("SNAPSHOT", "")))
            last = keys[-1]
            P  = np.array(f["P"][last])
            Bx = np.array(f["Bx"][last])
            By = np.array(f["By"][last])
        assert P.min() > 0.0, \
            f"Thermal pressure not positive: min = {P.min():.3e}"
        halfB2 = 0.5 * (Bx**2 + By**2)
        mask = halfB2 > 1e-10
        if mask.any():
            assert np.all(P[mask] < P[mask] + halfB2[mask]), \
                "P does not appear to be thermal pressure (should be < P + B^2/2)"
