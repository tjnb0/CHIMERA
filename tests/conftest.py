"""
Shared pytest fixtures for CHIMERA tests.

Fixtures:
  built_executable  -- compiles once per session, returns Path to chimera
  run_sim           -- function-scoped: run the binary, return HDF5 path
  ot_128_h5         -- session-scoped: run OT at N=128, return HDF5 path
                       (used by validation tests; only executes when needed)
"""

import os
import subprocess
import pytest
from pathlib import Path

REPO_ROOT  = Path(__file__).resolve().parent.parent
EXECUTABLE = REPO_ROOT / ("chimera.exe" if os.name == "nt" else "chimera")


def pytest_configure(config):
    """Register custom markers so pytest does not emit warnings."""
    config.addinivalue_line(
        "markers",
        "slow: marks tests as slow - OT vortex at N=128 "
        "(deselect with '-m \"not slow\"')"
    )


@pytest.fixture(scope="session")
def built_executable():
    """Build the main CHIMERA executable once for the entire session."""
    result = subprocess.run(
        ["make", "all"],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        pytest.fail(
            f"make failed - fix compilation errors before running tests.\n"
            f"--- stderr ---\n{result.stderr}"
        )
    if not EXECUTABLE.is_file():
        pytest.fail(f"Executable not found after make: {EXECUTABLE}")
    return EXECUTABLE


@pytest.fixture
def run_sim(built_executable, tmp_path):
    """
    Run the CHIMERA binary with given arguments, return path to HDF5 output.

    Usage:
        h5 = run_sim(problem_type=1, N=32)
        h5 = run_sim(problem_type=5, N=64, seed=42, h5_name="mc.h5")

        # GRF stress-test overrides (args 6 & 7 to the binary):
        #   extra_args=[target_M_s, target_beta]
        #   Pass -1 for either to keep random sampling for that parameter.
        h5 = run_sim(problem_type=5, N=64, seed=1,
                     extra_args=[5.0, -1.0])   # high-Mach, beta random
        h5 = run_sim(problem_type=5, N=64, seed=1,
                     extra_args=[-1.0, 0.1])   # M_s random, low-beta
    """
    def _run(problem_type, N, seed=0, h5_name=None, timeout=180,
             extra_args=None):
        if h5_name is None:
            h5_name = f"test_p{problem_type}_N{N}_s{seed}.h5"
        out_dir = str(tmp_path) + os.sep
        cmd = [
            str(built_executable), str(problem_type), str(N), str(seed),
            out_dir, h5_name,
        ]
        if extra_args is not None:
            cmd.extend(str(a) for a in extra_args)

        result = subprocess.run(
            cmd,
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if result.returncode != 0:
            pytest.fail(
                f"Simulation failed (problem={problem_type}, N={N}, seed={seed}, "
                f"extra_args={extra_args})\n"
                f"--- stdout ---\n{result.stdout}\n"
                f"--- stderr ---\n{result.stderr}"
            )
        h5_path = tmp_path / h5_name
        if not h5_path.exists():
            pytest.fail(f"HDF5 not found after simulation: {h5_path}\n"
                        f"stdout: {result.stdout}")
        return h5_path

    return _run


@pytest.fixture(scope="session")
def ot_128_h5(built_executable, tmp_path_factory):
    """
    Run OT vortex at N=128 exactly once per test session.
    Used by the validation tests; only executes if those tests are collected.
    """
    tmp     = tmp_path_factory.mktemp("validation")
    h5_name = "ot_N128.h5"
    result  = subprocess.run(
        [str(built_executable), "1", "128", "0",
         str(tmp) + os.sep, h5_name],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=600,
    )
    if result.returncode != 0:
        pytest.fail(
            f"OT N=128 simulation failed:\n"
            f"--- stdout ---\n{result.stdout}\n"
            f"--- stderr ---\n{result.stderr}"
        )
    return tmp / h5_name
