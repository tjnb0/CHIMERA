import subprocess
from pathlib import Path
from tqdm import tqdm
import os

# ---------------------------------------------------------------------------
# Repository layout (relative to this script in CHIMERA/scripts/)
#
#   CHIMERA/
#   ├── src/          – Fortran source
#   ├── scripts/      – this file
#   ├── obj/          – compiler objects   (auto-created)
#   ├── mod/          – compiled modules   (auto-created)
#   ├── outputs/      – simulation output  (auto-created)
#   ├── Makefile
#   └── mhd_sim[.exe]
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent        # CHIMERA/scripts/
REPO_ROOT  = SCRIPT_DIR.parent                      # CHIMERA/
OUTPUT_DIR = REPO_ROOT / "outputs"
EXECUTABLE = REPO_ROOT / ("mhd_sim.exe" if os.name == "nt" else "mhd_sim")

# Auto-create build and output directories so a fresh clone never crashes
for _dir in [OUTPUT_DIR, REPO_ROOT / "obj", REPO_ROOT / "mod"]:
    _dir.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------
N_RUNS       = 1000
SEED_START   = 1
TEST_PROBLEM = "5"


def run_one(run_id: int) -> None:
    """Execute the Fortran binary for a single Monte Carlo realisation."""
    seed    = SEED_START + run_id
    h5_name = f"mc_run_{seed}.h5"

    if not EXECUTABLE.is_file():
        raise FileNotFoundError(
            f"Executable not found: {EXECUTABLE}\n"
            f"Run 'make' in {REPO_ROOT} first."
        )

    subprocess.run(
        [str(EXECUTABLE), TEST_PROBLEM, str(seed),
         str(OUTPUT_DIR) + os.sep, h5_name],
        cwd=str(REPO_ROOT),   # run from repo root so Fortran relative paths are stable
        check=True,
        timeout=60,
    )


if __name__ == "__main__":
    print(f"\nRunning {N_RUNS} Monte Carlo simulations (seeds {SEED_START}–{SEED_START + N_RUNS - 1})")
    print(f"Output directory: {OUTPUT_DIR}\n")

    failed = 0
    with tqdm(total=N_RUNS) as pbar:
        for i in range(N_RUNS):
            try:
                run_one(i)
            except subprocess.CalledProcessError as e:
                print(f"\nRun {i + 1} failed (return code {e.returncode})")
                failed += 1
            except subprocess.TimeoutExpired:
                print(f"\nRun {i + 1} timed out")
                failed += 1
            finally:
                pbar.update(1)

    total = N_RUNS
    print(f"\nDone. {total - failed}/{total} runs completed successfully.")
    if failed:
        print(f"       {failed} run(s) failed - check output above.")
