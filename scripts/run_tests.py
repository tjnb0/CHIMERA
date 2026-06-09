#!/usr/bin/env python3
"""
CHIMERA unified test runner.

Usage:
    python scripts/run_tests.py                  # structural + physics + stress (default)
    python scripts/run_tests.py --structural     # Fortran unit tests + structural Python
    python scripts/run_tests.py --physics        # physics fidelity tests
    python scripts/run_tests.py --stress         # stress / GRF tests
    python scripts/run_tests.py --validation     # OT validation tests (slow, N=128)
    python scripts/run_tests.py --all            # everything

Exits 0 if all selected tests pass, 1 if any fail.

Validation tests require tests/validation/ot_reference.h5 to exist for the
L2 density and pressure profile checks. Generate it from an external trusted
source (Athena++, PLUTO, or published benchmark data) in the same CHIMERA
HDF5 format - one SNAPSHOT1 group per field at t=0.5, N=128.
"""

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT  = Path(__file__).resolve().parent.parent
TESTS_DIR  = REPO_ROOT / "tests"
FORT_DIR   = TESTS_DIR / "fortran"
FORT_TESTS = ["test_change_states", "test_field_ops", "test_bcs"]

PYTEST_DIRS = {
    "structural": TESTS_DIR / "structural",
    "physics":    TESTS_DIR / "physics",
    "stress":     TESTS_DIR / "stress",
    "validation": TESTS_DIR / "validation",
}


def run(cmd, cwd=None):
    return subprocess.run(cmd, cwd=str(cwd or REPO_ROOT)).returncode


def header(text):
    print(f"\n{'=' * 60}\n  {text}\n{'=' * 60}")


def build(include_fortran_tests):
    header("Building")
    targets = ["all"] + (["tests"] if include_fortran_tests else [])
    ret = run(["make"] + targets)
    if ret != 0:
        sys.exit("ERROR: Build failed. Fix compilation errors first.")
    print("Build OK.")


def run_fortran_tests():
    header("Fortran unit tests")
    n_pass = n_fail = n_skip = 0
    for name in FORT_TESTS:
        exe = FORT_DIR / name
        print(f"\n--- {name} ---")
        if not exe.is_file():
            print(f"  SKIP: not built ({exe})")
            n_skip += 1
            continue
        ret = run([str(exe)])
        n_pass += (ret == 0)
        n_fail += (ret != 0)
    print(f"\nFortran: {n_pass} passed, {n_fail} failed, {n_skip} skipped")
    return n_fail == 0


def run_pytest(dirs, slow=False):
    check = subprocess.run(
        [sys.executable, "-m", "pytest", "--version"],
        capture_output=True, text=True
    )
    if check.returncode != 0:
        print("ERROR: pytest not found.  pip install pytest")
        return False

    cmd = [sys.executable, "-m", "pytest", "--tb=short", "--no-header", "-v"]
    cmd += [str(d) for d in dirs if d.exists()]
    if not slow:
        cmd += ["-m", "not slow"]
    return run(cmd) == 0


def main():
    parser = argparse.ArgumentParser(description="CHIMERA test runner")
    parser.add_argument("--structural", action="store_true")
    parser.add_argument("--physics",    action="store_true")
    parser.add_argument("--stress",     action="store_true")
    parser.add_argument("--validation", action="store_true",
                        help="Include slow OT N=128 validation tests")
    parser.add_argument("--all",        action="store_true",
                        help="Run everything including validation")
    args = parser.parse_args()

    run_default = not any([args.structural, args.physics,
                           args.stress, args.validation, args.all])

    do_structural = args.all or args.structural or run_default
    do_physics    = args.all or args.physics    or run_default
    do_stress     = args.all or args.stress     or run_default
    do_validation = args.all or args.validation

    build(include_fortran_tests=do_structural)

    results = {}

    if do_structural:
        header("Structural tests")
        fort_ok   = run_fortran_tests()
        pytest_ok = run_pytest([PYTEST_DIRS["structural"]])
        results["structural"] = fort_ok and pytest_ok

    if do_physics:
        header("Physics fidelity tests")
        results["physics"] = run_pytest([PYTEST_DIRS["physics"]])

    if do_stress:
        header("Stress tests")
        results["stress"] = run_pytest([PYTEST_DIRS["stress"]])

    if do_validation:
        header("Validation tests (slow)")
        results["validation"] = run_pytest(
            [PYTEST_DIRS["validation"]], slow=True
        )

    header("Summary")
    all_pass = True
    for name, ok in results.items():
        status = "PASS" if ok else "FAIL"
        print(f"  {name:<16} {status}")
        all_pass = all_pass and ok

    print()
    if all_pass:
        print("All tests passed.")
        sys.exit(0)
    else:
        print("Some tests failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
