import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from datetime import datetime
from pathlib import Path
import re
import h5py

try:
    #--------------------
    #--- Main Options ---
    #--------------------
    test_problem = "5"   # 1 = Orszag-Tang, 2 = KH, 3 = Field Loop, 4 = Rotor, 5 = MC
    fortran_seed = 42   # 0 = primitive_snaps_MC.h5, else use desired MC run number
    run_new_sim  = True  # Create new data or plot existing data


    #----------------------
    # --- Other Options ---
    #----------------------
    save_gif    = False         # Save time evolution as GIF
    save_plots  = False         # Save final frame as PNG
    output_path  = "../outputs/"

    if (fortran_seed == 0):
        h5_filename = "primitive_snaps_MC.h5"
    else:
        h5_filename = "primitive_snaps_MC"+str(fortran_seed)+".h5"


    #-----------------------------
    #--- Get Test Problem Name ---
    #-----------------------------
    title, pic_name = {
        "1": ("Orszag-Tang Vortex: "          , "Orszag_Tang_Vortex_"),
        "2": ("Kelvin-Helmholtz Instability: ", "Kelvin_Helmholtz_Instability_"),
        "3": ("Field Loop Advection: "        , "Field_Loop_Advection_"),
        "4": ("Periodic MHD Rotor: "          , "Periodic_MHD_Rotor_"),
        "5": ("Gaussian Random Field: "       , "Gaussian_Random_Field"),
    }.get(test_problem, (None, None))
    if title is None:
        print("Undefined test_problem:", test_problem)
        exit()


    #-------------------------
    # --- Helper Functions ---
    #-------------------------
    def read_N_from_config(filename="src/mhd_config.f90"):
        with open(filename, "r") as f:
            for line in f:
                line = line.split('!')[0]
                m = re.match(r"\s*integer\s*,\s*parameter\s*::\s*N\s*=\s*(\d+)",
                             line, re.IGNORECASE)
                if m:
                    return int(m.group(1))
        raise ValueError("N not found in config file.")


    def load_hdf5_snapshots(h5_path):
        """
        Open the HDF5 file and return all 6 fields as stacked snapshot arrays.
        Layout: /rho/SNAPSHOT1 ..., /P/..., /Bx/..., /By/..., /Vx/..., /Vy/...
                /time/time_array
        """
        f = h5py.File(h5_path, "r")
        data = {}
        for field in ["rho", "P", "Bx", "By", "Vx", "Vy"]:
            if field in f:
                group = f[field]
                snaps = [group[k][:] for k in sorted(group.keys(),
                         key=lambda x: int(x.replace("SNAPSHOT", "")))]
                data[field] = np.stack(snaps)
            else:
                data[field] = None
        data["time"] = f["time/sim_time"][:] if "time" in f else None
        data["file"] = f
        return data


    def plot_all_fields(h5_path, suptitle, output_image, save_fig, save_gif):
        """
        Plot all 6 MHD fields on a single 2x3 figure, animating through snapshots.
        Layout:
            row 0: rho  |  Bx  |  By
            row 1: P    |  Vx  |  Vy
        """
        print("Plotting all fields...")

        data = load_hdf5_snapshots(h5_path)
        time_arr = data.get("time", None)

        # Field layout: (key, label)
        layout = [
            [("rho", r"$\rho$"),  ("Bx", r"$B_x$"),  ("By", r"$B_y$")],
            [("P",   r"$P$"),     ("Vx", r"$V_x$"),  ("Vy", r"$V_y$")],
        ]

        # Precompute global color limits for each field
        global_vmin = {}
        global_vmax = {}
        for row in range(2):
            for key, _ in layout[row]:
                arr = data[key]          # shape: (n_snap, Ny, Nx)
                global_vmin[key] = arr.min()
                global_vmax[key] = arr.max()

        n_snap = data["rho"].shape[0]
        extent = [0, 1, 0, 1]
        ticks  = [0, 0.25, 0.5, 0.75, 1]
        labels = ['0', '0.25', '0.5', '0.75', '1']

        fig, axes = plt.subplots(2, 3, figsize=(13, 8))
        fig.suptitle(suptitle, fontsize=13)
        ims = []; keys = []; cbs = []

        # --- Initialise with first snapshot ---
        for row in range(2):
            for col in range(3):
                key, lbl = layout[row][col]
                ax = axes[row, col]
                frame = data[key][0]

                im = ax.imshow(frame, origin="lower", aspect="auto",
                            cmap="viridis", extent=extent,
                            vmin=global_vmin[key], vmax=global_vmax[key])

                cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

                ax.set_title(lbl, fontsize=11)
                ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=7)
                ax.set_yticks(ticks); ax.set_yticklabels(labels, fontsize=7)

                ims.append(im)
                keys.append(key)
                cbs.append(cb)

        t_str = f"t = {time_arr[0]:.2f}" if time_arr is not None else "t = ?"
        time_text = fig.text(0.5, 0.96, t_str, ha="center", va="top", fontsize=11)

        # --- Animate interactively ---
        plt.ion()
        try:
            def update(i):
                artists = []

                for im, key in zip(ims, keys):
                    im.set_data(data[key][i])
                    artists.append(im)

                t_str = f"t = {time_arr[i]:.3f}" if time_arr is not None else f"snap {i}"
                time_text.set_text(t_str)
                artists.append(time_text)

                return artists

            ani = animation.FuncAnimation(
                fig,
                update,
                frames=n_snap,
                interval=1,
                blit=True
            )
            plt.show()

        except KeyboardInterrupt:
            print("\nAnimation interrupted by user.")

        plt.ioff()

        # --- Save final frame as PNG ---
        if test_problem == "5":
            if fortran_seed == 0:
                suffix = "_MC"
            else:
                suffix = f"_MC{fortran_seed}"
        else:
            suffix = ""
        if save_fig:
            # Force final frame before saving
            final_i = n_snap - 1
            for im, key in zip(ims, keys):
                im.set_data(data[key][final_i])
            if time_arr is not None:
                time_text.set_text(f"t = {time_arr[final_i]:.3f}")
            else:
                time_text.set_text(f"snap {final_i}")
            fig.canvas.draw()
            png_path = f"{output_image}{suffix}.png"
            fig.savefig(png_path, dpi=240)
            print(f"Saved PNG: {png_path}")

        # --- Save as GIF ---
        if save_gif:
            def update(i):
                artists = []
                for im, key in zip(ims, keys):
                    frame = data[key][i]
                    im.set_data(frame)
                    artists.append(im)
                t_str = f"t = {time_arr[i]:.3f}" if time_arr is not None else f"snap {i}"
                time_text.set_text(t_str)
                artists.append(time_text)
                return artists

            ani = animation.FuncAnimation(fig, update, frames=n_snap,
                                          blit=False, repeat=False)
            gif_path = f"{output_image}{suffix}.gif"
            ani.save(gif_path, writer=animation.PillowWriter(fps=10))
            print(f"Saved GIF: {gif_path}")

        plt.show()
        data["file"].close()


    #-----------------------------
    #--- Run Fortran Sim: T/F? ---
    #-----------------------------
    if run_new_sim:

        output_dir = "outputs"
        os.makedirs(output_dir, exist_ok=True)
        print(f"\nDeleting old *.h5 files in ./{output_dir}/")
        for filename in os.listdir(output_dir):
            if filename.endswith(".h5"):
                os.remove(os.path.join(output_dir, filename))

        print("\nStarting " + title[:-2] + " simulation.")
        start_time = datetime.now()
        print("Start Time:", start_time.strftime("%I:%M:%S %p %Z %Y-%m-%d"), "\n")

        if os.name == "nt":
            executable = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mhd_sim.exe")
        else:
            executable = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mhd_sim")

        if not os.path.isfile(executable):
            print(f"Error: The executable '{executable}' does not exist.")
        else:
            subprocess.run([executable, test_problem, str(fortran_seed), output_path, h5_filename])

        print("Simulation complete.")
        end_time = datetime.now()
        print("End Time:", end_time.strftime("%I:%M:%S %p %Z %Y-%m-%d"))
        print("Elapsed :", str(end_time - start_time), "\n")

    else:
        print("\nReading old data...")


    #----------------------
    # --- Main Function ---
    #----------------------
    def main():
        N = read_N_from_config()
        output_dir = Path(__file__).parent / "outputs/"
        h5_path    = output_dir / h5_filename

        plot_all_fields(
            h5_path,
            suptitle     = title.strip(": "),
            output_image = output_dir / pic_name,
            save_fig     = save_plots,
            save_gif     = save_gif,
        )

        print("Exiting program.\n")

    if __name__ == "__main__":
        main()


except KeyboardInterrupt:
    print("\nExecution interrupted by user (Ctrl+C).")
