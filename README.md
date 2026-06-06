# CHIMERA
### Compressible, High-resolution Ideal MHD with Ensemble Randomisation Approach

A finite-volume solver for the compressible ideal magnetohydrodynamics (MHD)
equations in two spatial dimensions, written in modern Fortran. Developed as
an independent project.

---

## Overview

CHIMERA evolves the compressible ideal MHD equations in conservation form on
a uniform 2D Cartesian grid. The numerical scheme combines second-order
MUSCL-Hancock reconstruction with a Rusanov (local Lax-Friedrichs) Riemann
solver, and uses constrained transport to preserve the divergence-free
condition on B to machine precision. OpenMP threading accelerates the
reconstruction and slope-limiting passes. Boundary conditions are configurable
per side, supporting periodic, zero-gradient outflow, fixed, and driven inflow
on each of the four domain edges independently.

Output is written to HDF5, with each field stored as a sequence of snapshots
alongside realised physics diagnostics (gamma, Mach number, plasma beta).

---

## Governing Equations

CHIMERA advances the compressible MHD system in conservation form:

```
d(rho)/dt   + div(rho*v)                       = 0      (mass)
d(rho*v)/dt + div(rho*v*v + P*I - B*B)         = 0      (momentum)
d(E)/dt     + div((E + P*)*v - B*(v.B))        = 0      (energy)
d(B)/dt     - curl(v x B)                      = 0      (induction)
```

where `P* = p + |B|^2/2` is the total (thermal + magnetic) pressure and `E`
is the total energy density (internal + kinetic + magnetic). The gas is closed
by an ideal equation of state with adiabatic index gamma.

---

## Numerical Methods

| Component | Method |
|-----------|--------|
| Spatial discretisation | Cell-centred finite volume on a uniform Cartesian grid |
| Time integration | Predictor-corrector (MUSCL-Hancock); CFL-limited adaptive timestep |
| Reconstruction | 2nd-order MUSCL with MOOD fallback to 1st-order at troubled cells |
| Slope limiting | Van Leer harmonic mean limiter (toggleable) |
| Riemann solver | Local Lax-Friedrichs / Rusanov |
| Divergence control | Constrained transport (CT) on staggered face-centred B; div B monitored every step |
| Parallelism | OpenMP on reconstruction and slope-limiting loops |
| Boundary conditions | Per-side ghost-cell layer: periodic, outflow, fixed, or driven inflow |

### Time-stepping

At each step CHIMERA:
1. Fills ghost cells for all six primitive fields using the per-side BC flags.
2. Computes BC-aware gradients and applies the Van Leer slope limiter.
3. Predicts primitive variables half a timestep forward (MUSCL-Hancock prediction).
4. Reconstructs left/right face states with MOOD fallback and a thermal pressure
   positivity check to prevent unphysical states in high-field regions.
5. Evaluates Rusanov fluxes and updates conserved variables.
6. Advances face-centred B via constrained transport.

---

## Features

- **MUSCL-Hancock predictor-corrector** - second-order accurate in space and time
- **MOOD reconstruction** - per-cell fallback to first order where reconstructed
  values exceed stencil bounds or implied thermal pressure falls below the floor
- **Van Leer slope limiter** - toggled via `useSlopeLimiting` in `mhd_config.f90`
- **Constrained transport** - staggered face-centred B updated via discrete curl
  of Ez; div B monitored and printed each timestep
- **Fast magnetosonic CFL condition** - timestep limited by `c_f + |v|` with a
  Courant factor of 0.3; correct thermal pressure used for sound speed
- **Generic open boundary conditions** - four independent per-side flags
  (`BC_xlo`, `BC_xhi`, `BC_ylo`, `BC_yhi`) with ghost-cell padding each step;
  supports periodic, zero-gradient outflow, fixed, and driven inflow (CME-ready)
- **Gaussian random field (GRF) initial conditions** - divergence-free velocity
  (from stream function) and magnetic field (from vector potential) generated
  from RBF-smoothed random fields; scaled to target V0 and B0 with a fast-speed
  cap; realised gamma, sonic Mach number, and plasma beta stored in output
- **HDF5 output** - all primitive fields written as snapshot groups with a shared
  time array and physics parameters; pressure output is thermal pressure p
- **OpenMP threading** - slope-limiting and MOOD reconstruction loops
  parallelised with `!$omp parallel do`

---

## Test Problems

CHIMERA ships with five built-in initial conditions selected via command-line:

| `problem_type` | Problem | BCs | Physics tested |
|:-:|---------|-----|----------------|
| 1 | **Orszag-Tang vortex** | Periodic | 2D MHD turbulence, shock-current-sheet interaction |
| 2 | **Kelvin-Helmholtz instability** | Periodic | Shear-driven magnetic instability and interface roll-up |
| 3 | **Field-loop advection** | Periodic | Accuracy of CT scheme; passive advection of a magnetic flux tube |
| 4 | **MHD rotor** | Outflow (all sides) | High-density rotating disk; torsional Alfven waves, open boundaries |
| 5 | **Monte Carlo / GRF ensemble** | Periodic | Statistical studies of compressible MHD turbulence with randomised ICs |

Recommended resolutions: Orszag-Tang N=3000, KH N=1024, Advection N=1024,
Rotor N=2000.

---

## Code Structure

```
.
|-- src/
|   |-- main.f90                 - Time loop, CFL timestep, I/O scheduling,
|   |                              outflow boundary flux correction
|   |-- mhd_config.f90           - Global parameters (N, tEnd, CFL, floors,
|   |                              BC type constants and per-side flags)
|   |-- mhd_init.f90             - Array allocation, grid setup, all five
|   |                              initial conditions including GRF generation
|   |-- mhd_bc.f90               - Ghost-cell BC module: fills (N+2)x(N+2)
|   |                              padded arrays for each primitive field
|   |-- mhd_change_states.f90    - Primitive <-> conserved variable conversion
|   |-- mhd_derivatives.f90      - BC-aware gradients, Van Leer slope limiter,
|   |                              MUSCL/MOOD face reconstruction,
|   |                              thermal pressure positivity check
|   |-- mhd_flux.f90             - Rusanov flux evaluation, BC-aware conserved
|   |                              variable update, constrained transport
|   |-- mhd_field_ops.f90        - BC-aware discrete curl, div B diagnostic,
|   |                              face-to-cell B averaging
|   |-- mhd_write_h5.f90         - HDF5 output (snapshots + physics parameters)
|-- scripts/
|   |-- run_mhd.py               - Run single simulation and animate output
|   `-- run_mhd_MC.py            - Run Monte Carlo ensemble
|-- Makefile
`-- README.md
```

---

## Getting Started

### Prerequisites

- Fortran compiler: `gfortran >= 9` or Intel `ifx`/`ifort`
- HDF5 library with Fortran bindings (e.g. `libhdf5-fortran-dev` on Debian/Ubuntu)
- OpenMP (included with most compilers)
- Python 3 with `h5py`, `numpy`, `matplotlib` (for scripts)

### Build

```bash
git clone https://github.com/tjnb0/CHIMERA.git
cd CHIMERA
make
```

Edit the `Makefile` to point to your HDF5 installation if needed.

### Run

```bash
# Single run - Orszag-Tang vortex
./mhd_sim 1 0 ./outputs/ orszag_tang.h5

# MHD Rotor with outflow BCs
./mhd_sim 4 0 ./outputs/ rotor.h5

# Monte Carlo ensemble (100 runs)
for seed in $(seq 1 100); do
    ./mhd_sim 5 $seed ./outputs/ mc_run_${seed}.h5
done
```

Or using the Python scripts from the `scripts/` directory:

```bash
python scripts/run_mhd.py          # single run with interactive animation
python scripts/run_mhd_MC.py       # full Monte Carlo ensemble
```

**Command-line arguments:** `problem_type  seed  output_path  filename.h5`

Grid resolution, end time, and scheme options are compile-time parameters
in `mhd_config.f90`.

---

## Boundary Conditions

Each side of the domain is assigned independently in the problem setup routine:

```fortran
BC_xlo = BC_PERIODIC   ! left
BC_xhi = BC_PERIODIC   ! right
BC_ylo = BC_PERIODIC   ! bottom
BC_yhi = BC_PERIODIC   ! top
```

Available types defined in `mhd_config.f90`:

| Constant | Value | Behaviour |
|----------|-------|-----------|
| `BC_PERIODIC` | 1 | Circular wrap (default for all periodic problems) |
| `BC_OUTFLOW` | 2 | Zero-gradient outflow; wave exits cleanly |
| `BC_FIXED` | 3 | Prescribed ambient state via `f_ambient` argument |
| `BC_INFLOW` | 4 | Driven inflow; intended for CME-type problems |

---

## Output Format

Results are written to a single HDF5 file per run:

```
/rho/SNAPSHOT1 ... SNAPSHOTn     - density              [N x N, float64]
/P/SNAPSHOT1   ... SNAPSHOTn     - thermal pressure p    [N x N, float64]
/Bx/SNAPSHOT1  ... SNAPSHOTn     - x magnetic field
/By/SNAPSHOT1  ... SNAPSHOTn     - y magnetic field
/Vx/SNAPSHOT1  ... SNAPSHOTn     - x velocity
/Vy/SNAPSHOT1  ... SNAPSHOTn     - y velocity
/time/sim_time                    - snapshot times        [n, float64]
/parameters/gamma                 - adiabatic index (realised)
/parameters/M_s                   - sonic Mach number (realised)
/parameters/beta                  - plasma beta = 2p/B^2 (realised)
```

Note: the `/P/` dataset contains **thermal pressure** `p = P* - B^2/2`, not
total pressure. Snapshots are written every `tOut = 0.01` time units.

```python
import h5py, numpy as np
with h5py.File("orszag_tang.h5", "r") as f:
    rho  = np.array(f["rho/SNAPSHOT10"])
    p    = np.array(f["P/SNAPSHOT10"])     # thermal pressure
    time = np.array(f["time/sim_time"])
```

---

## References

- Toth, G. (2000). *The div(B)=0 constraint in shock-capturing MHD codes.* J. Comput. Phys.
- Evans, C. R. & Hawley, J. F. (1988). *Simulation of magnetohydrodynamic flows:
  A constrained transport method.* ApJ.
- van Leer, B. (1979). *Towards the ultimate conservative difference scheme.*
  J. Comput. Phys.
- Clain, S., Diot, S. & Loubere, R. (2011). *A high-order finite volume method
  for hyperbolic systems: Multi-dimensional Optimal Order Detection (MOOD).*
  J. Comput. Phys.

---

## License

MIT License - see [LICENSE](LICENSE) for details.
