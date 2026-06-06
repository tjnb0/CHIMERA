import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from datetime import datetime
from pathlib import Path
import h5py

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

SCRIPT_DIR  = Path(__file__).resolve().parent       # CHIMERA/scripts/
REPO_ROOT   = SCRIPT_DIR.parent                     # CHIMERA/
SRC_DIR     = REPO_ROOT / "src"
OUTPUT_DIR  = REPO_ROOT / "outputs"
EXECUTABLE  = REPO_ROOT / ("mhd_sim.exe" if os.name == "nt" else "mhd_sim")

# Auto-create build and output directories so a fresh clone never crashes
for _dir in [OUTPUT_DIR, REPO_ROOT / "obj", REPO_ROOT / "mod"]:
    _dir.mkdir(parents=True, exist_ok=True)


try:
    # -----------------------------------------------------------------------
    # Main options
    # -----------------------------------------------------------------------
    test_problem = "5"   # 1=Orszag-Tang, 2=KH, 3=Field Loop, 4=Rotor, 5=MC
    N            = 64    # Grid size (N x N cells)
    fortran_seed = 42    # Seed for MC runs (ignored for problems 1-4)
    run_new_sim  = True  # True=run simulation, False=plot existing output

    # -----------------------------------------------------------------------
    # Output options
    # -----------------------------------------------------------------------
    save_gif   = False
    save_plots = False

    # -----------------------------------------------------------------------
    # Derive filenames from problem type
    # -----------------------------------------------------------------------
    PROBLEM_META = {
        "1": ("Orszag-Tang Vortex",           "Orszag_Tang_Vortex",          "orszag_tang.h5"),
        "2": ("Kelvin-Helmholtz Instability", "Kelvin_Helmholtz_Instability","kelvin_helmholtz.h5"),
        "3": ("Field Loop Advection",         "Field_Loop_Advection",        "field_loop.h5"),
        "4": ("Periodic MHD Rotor",           "Periodic_MHD_Rotor",          "mhd_rotor.h5"),
        "5": ("Gaussian Random Field",        "Gaussian_Random_Field",       None),  # set below
    }

    if test_problem not in PROBLEM_META:
        print(f"Undefined test_problem: {test_problem}")
        raise SystemExit

    title_str, pic_name, h5_filename = PROBLEM_META[test_problem]

    if test_problem == "5":
        h5_filename = "mc_run.h5" if fortran_seed == 0 else f"mc_run_{fortran_seed}.h5"

    h5_path = OUTPUT_DIR / h5_filename

    # -----------------------------------------------------------------------
    # Helper: load all primitive-variable snapshots from HDF5
    # -----------------------------------------------------------------------
    def load_hdf5_snapshots(path: Path) -> dict:
        f = h5py.File(path, "r")
        data = {}
        for field in ["rho", "P", "Bx", "By", "Vx", "Vy"]:
            if field in f:
                group = f[field]
                snaps = [
                    group[k][:]
                    for k in sorted(group.keys(),
                                    key=lambda x: int(x.replace("SNAPSHOT", "")))
                ]
                data[field] = np.stack(snaps)
            else:
                data[field] = None
        data["time"] = f["time/sim_time"][:] if "time" in f else None
        data["file"] = f
        return data

    # -----------------------------------------------------------------------
    # Helper: animate all 6 MHD fields; optionally save PNG / GIF
    # -----------------------------------------------------------------------
    def plot_all_fields(path: Path, suptitle: str, output_stem: Path,
                        save_fig: bool, save_gif: bool) -> None:
        print("Plotting all fields...")
        data = load_hdf5_snapshots(path)
        time_arr = data.get("time")

        layout = [
            [("rho", r"$\rho$"), ("Bx", r"$B_x$"), ("By", r"$B_y$")],
            [("P",   r"$P$"),    ("Vx", r"$V_x$"), ("Vy", r"$V_y$")],
        ]

        # Global colour limits per field
        vmin = {key: data[key].min() for row in layout for key, _ in row}
        vmax = {key: data[key].max() for row in layout for key, _ in row}

        n_snap = data["rho"].shape[0]
        extent = [0, 1, 0, 1]
        ticks  = [0, 0.25, 0.5, 0.75, 1.0]
        labels = ["0", "0.25", "0.5", "0.75", "1"]

        fig, axes = plt.subplots(2, 3, figsize=(13, 8))
        fig.suptitle(suptitle, fontsize=13)
        ims = []; keys_list = []

        for row in range(2):
            for col in range(3):
                key, lbl = layout[row][col]
                ax = axes[row, col]
                im = ax.imshow(
                    data[key][0], origin="lower", aspect="auto",
                    cmap="viridis", extent=extent,
                    vmin=vmin[key], vmax=vmax[key]
                )
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                ax.set_title(lbl, fontsize=11)
                ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=7)
                ax.set_yticks(ticks); ax.set_yticklabels(labels, fontsize=7)
                ims.append(im); keys_list.append(key)

        t_str     = f"t = {time_arr[0]:.2f}" if time_arr is not None else "t = ?"
        time_text = fig.text(0.5, 0.96, t_str, ha="center", va="top", fontsize=11)

        def _update(i):
            for im, key in zip(ims, keys_list):
                im.set_data(data[key][i])
            t_label = f"t = {time_arr[i]:.3f}" if time_arr is not None else f"snap {i}"
            time_text.set_text(t_label)
            return ims + [time_text]

        # --- Interactive animation ---
        plt.ion()
        try:
            # blit=False: fig.text() is a figure-level artist (no parent axes),
            # so blit=True causes a Tkinter/matplotlib crash on ax._get_view().
            # interval=50 → smooth 20 fps; interval=1 was an unreachable 1000 fps.
            ani = animation.FuncAnimation(        # noqa: F841
                fig, _update, frames=n_snap, interval=50, blit=False
            )
            plt.show()
        except KeyboardInterrupt:
            print("\nAnimation interrupted by user.")
        plt.ioff()

        # --- Save final frame as PNG ---
        if save_fig:
            final = n_snap - 1
            _update(final)
            fig.canvas.draw()
            png_path = output_stem.with_suffix(".png")
            fig.savefig(png_path, dpi=240)
            print(f"Saved PNG: {png_path}")

        # --- Save as GIF ---
        if save_gif:
            ani_gif = animation.FuncAnimation(
                fig, _update, frames=n_snap, blit=False, repeat=False
            )
            gif_path = output_stem.with_suffix(".gif")
            ani_gif.save(gif_path, writer=animation.PillowWriter(fps=10))
            print(f"Saved GIF: {gif_path}")

        plt.show()
        data["file"].close()

    # -----------------------------------------------------------------------
    # Optionally run the Fortran simulation
    # -----------------------------------------------------------------------
    if run_new_sim:
        print(f"\nClearing old *.h5 files in {OUTPUT_DIR}")
        for f in OUTPUT_DIR.glob("*.h5"):
            f.unlink()

        print(f"\nStarting {title_str} simulation.")
        start = datetime.now()
        print("Start:", start.strftime("%I:%M:%S %p  %Y-%m-%d"), "\n")

        if not EXECUTABLE.is_file():
            raise FileNotFoundError(
                f"Executable not found: {EXECUTABLE}\n"
                f"Run 'make' in {REPO_ROOT} first."
            )

        subprocess.run(
            [str(EXECUTABLE), test_problem, str(N), str(fortran_seed),
             str(OUTPUT_DIR) + os.sep, h5_filename],
            cwd=str(REPO_ROOT),   # run from repo root so relative paths in Fortran are stable
            check=True
        )

        end = datetime.now()
        print("\nSimulation complete.")
        print("End    :", end.strftime("%I:%M:%S %p  %Y-%m-%d"))
        print("Elapsed:", str(end - start), "\n")

    else:
        print("\nReading existing data...")

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    def main() -> None:
        # Read N from the HDF5 output shape rather than the source file
        with h5py.File(h5_path, "r") as f:
            N_out = f["rho"][list(f["rho"].keys())[0]].shape[0]
        print(f"Grid size: N = {N_out}")

        plot_all_fields(
            path        = h5_path,
            suptitle    = title_str,
            output_stem = OUTPUT_DIR / pic_name,
            save_fig    = save_plots,
            save_gif    = save_gif,
        )
        print("Exiting.\n")

    if __name__ == "__main__":
        main()

except KeyboardInterrupt:
    print("\nExecution interrupted by user (Ctrl+C).")
