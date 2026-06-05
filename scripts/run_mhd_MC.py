import subprocess
from tqdm import tqdm
import os


# --------------
# --- Config ---
# --------------
N_RUNS      = 1000
SEED_START  = 1


# Run MHD sim for a given seed value
TEST_PROBLEM = "5"
OUTPUT_PATH  = "../outputs/"

def run_one(run_id: int) -> str:
    """
    Execute the Fortran binary for a single Monte-Carlo realization.
    Returns a short summary string.
    """
    seed = SEED_START + run_id            
    h5_name = f"primitive_snaps_MC{run_id + 1}.h5"

    exe_name = "mhd_sim.exe" if os.name == "nt" else "mhd_sim"
    executable = os.path.join(os.path.dirname(os.path.abspath(__file__)), exe_name)

    cmd = [
        executable,
        TEST_PROBLEM,
        str(seed),
        OUTPUT_PATH,
        h5_name
    ]

    subprocess.run(cmd, check=True, timeout=60)

    return f"Run {run_id + 1}: seed={seed}, file={h5_name}. Done."


if __name__ == "__main__":
    print(f'\n Running {N_RUNS} simulations. Starting seed: {SEED_START}')
    with tqdm(total=N_RUNS) as pbar:
        for i in range(N_RUNS):
            try:
                # msg = run_one(i)
                # print(msg)
                run_one(i)
            except subprocess.CalledProcessError as e:
                # Option to log errors if needed
                print(f"Run {i+1} failed (return code {e.returncode})")
                pass
            except subprocess.TimeoutExpired:
                print(f"Run {i+1} timed out")
                pass
            finally:
                pbar.update(1)
