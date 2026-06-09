# CHIMERA
### Constrained-transport High-resolution Ideal MHD for Ensemble Randomization and Analysis

A finite-volume solver for the compressible ideal magnetohydrodynamics (MHD)
equations in two spatial dimensions, written in modern Fortran. Developed as
an independent project.

---

## Overview

CHIMERA evolves the compressible ideal MHD equations in conservative form on
a uniform 2D Cartesian grid. The numerical scheme combines second-order
MUSCL-Hancock reconstruction with a Rusanov (local Lax-Friedrichs) Riemann
solver and a monotonized central (MC) slope limiter, and uses constrained
transport to preserve the divergence-free condition on B to machine precision.
The limiter is selectable at build time (MC or Van Leer) via a single constant
in `mhd_config.f90`. OpenMP threading accelerates the reconstruction and
slope-limiting passes. Boundary conditions are configurable per side,
supporting periodic, zero-gradient outflow, fixed, and driven inflow on each
of the four domain edges independently.

Output is written to HDF5, with each field stored as a sequence of snapshots
alongside the realized physics parameters (gamma, Mach number, plasma beta).

---

## Example Outputs
<img width="648" height="509" alt="Orszag_Tang_Vortex_Density" src="https://github.com/user-attachments/assets/9f13b17a-f63f-40b4-8e13-7e484bba0d00" />

---

## Governing Equations

CHIMERA advances the compressible MHD system as:

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
| Spatial discretization | Cell-centered finite volume on a uniform Cartesian grid |
| Time integration | Predictor-corrector (MUSCL-Hancock); CFL-limited adaptive timestep |
| Reconstruction | 2nd-order MUSCL with MOOD fallback to 1st-order at troubled cells |
| Slope limiting | Monotonized central (default) or Van Leer mean limiter |
| Riemann solver | Local Lax-Friedrichs / Rusanov |
| Divergence control | Constrained transport (CT) on staggered face-centered B; div B monitored every step |
| Parallelism | OpenMP on reconstruction and slope-limiting loops |
| Boundary conditions | Per-side ghost-cell layer: periodic, outflow, fixed, or driven inflow |

### Time-stepping

At each step CHIMERA:
1. Fills ghost cells for all six primitive fields using the per-side BC flags.
2. Computes BC-aware gradients and applies the slope limiter (MC or Van Leer,
   selected at build time in `mhd_config.f90`).
3. Predicts primitive variables half a timestep forward (MUSCL-Hancock prediction).
4. Reconstructs left/right face states with MOOD fallback and a thermal pressure
   positivity check to prevent unphysical states in high-field regions.
5. Evaluates Rusanov fluxes and updates conserved variables.
6. Advances face-centered B via constrained transport.

---

## Features

- **MUSCL-Hancock predictor-corrector** - second-order accurate in space and time
- **MOOD reconstruction** - per-cell fallback to first order where reconstructed
  values exceed stencil bounds or implied thermal pressure falls below the floor
- **Slope limiter** - monotonized central (MC, default) or Van Leer mean; selected
  at build time by setting `slope_limiter` in `mhd_config.f90` and recompiling.
  MC is less diffusive and fully TVD; Van Leer is more conservative near strong shocks.
- **Constrained transport** - staggered face-centered B updated via discrete curl
  of Ez; div B monitored and printed each timestep
- **Fast magnetosonic CFL condition** - timestep limited by `c_f + |v|` with a
  Courant factor of 0.3; correct thermal pressure used for sound speed
- **Generic open boundary conditions** - four independent per-side flags
  (`BC_xlo`, `BC_xhi`, `BC_ylo`, `BC_yhi`) with ghost-cell padding each step;
  supports periodic, zero-gradient outflow, fixed, and driven inflow (CME-ready)
- **Gaussian random field (GRF) initial conditions** - divergence-free velocity
  (from stream function) and magnetic field (from vector potential) generated
  from RBF-smoothed random fields; scaled to target V0 and B0 with a fast-speed
  cap; realized gamma, sonic Mach number, and plasma beta stored in output
- **HDF5 output** - all primitive fields written as snapshot groups with a shared
  time array and physics parameters; pressure output is thermal pressure p
- **OpenMP threading** - slope-limiting and MOOD reconstruction loops
  parallelized with `!$omp parallel do`

---

## Test Problems

CHIMERA has five built-in initial conditions selected via command-line:

| `problem_type` | Problem | BCs | Physics tested |
|:-:|---------|-----|----------------|
| 1 | **Orszag-Tang vortex** | Periodic | 2D MHD turbulence, shock-current-sheet interaction |
| 2 | **Kelvin-Helmholtz instability** | Periodic | Shear-driven magnetic instability and interface roll-up |
| 3 | **Field-loop advection** | Periodic | Accuracy of CT scheme; passive advection of a magnetic flux tube |
| 4 | **MHD rotor** | Outflow (all sides) | High-density rotating disk; torsional Alfven waves, open boundaries |
| 5 | **Monte Carlo / GRF ensemble** | Periodic | Statistical studies of compressible MHD turbulence with randomized ICs |

---

## Code Structure

```
.
|-- src/
|   |-- main.f90                 - Time loop, CFL timestep, I/O scheduling,
|   |                              outflow boundary flux correction
|   |-- mhd_config.f90           - Global parameters (N, tEnd, CFL, floors,
|   |                              slope limiter selection, BC type constants)
|   |-- mhd_init.f90             - Array allocation, grid setup, all five
|   |                              initial conditions including GRF generation
|   |-- mhd_bc.f90               - Ghost-cell BC module: fills (N+2)x(N+2)
|   |                              padded arrays for each primitive field
|   |-- mhd_change_states.f90    - Primitive <-> conserved variable conversion
|   |-- mhd_derivatives.f90      - BC-aware gradients, MC/Van Leer slope limiter,
|   |                              MUSCL/MOOD face reconstruction,
|   |                              thermal pressure positivity check
|   |-- mhd_flux.f90             - Rusanov flux evaluation, BC-aware conserved
|   |                              variable update, constrained transport
|   |-- mhd_field_ops.f90        - BC-aware discrete curl, div B diagnostic,
|   |                              face-to-cell B averaging
|   |-- mhd_write_h5.f90         - HDF5 output (snapshots + physics parameters)
|-- tests/
|   |-- conftest.py              - Shared pytest fixtures (built_executable,
|   |                              run_sim, ot_128_h5)
|   |-- fortran/                 - Fortran unit tests; built with 'make tests'
|   |   |-- test_bcs.f90         - Ghost-cell BC correctness
|   |   |-- test_change_states.f90 - Primitive <-> conserved roundtrip
|   |   `-- test_field_ops.f90   - Curl, div B, and gradient operators
|   |-- physics/                 - Physics fidelity tests
|   |   `-- test_physics.py      - Mass conservation, positivity, field loop
|   |                              convergence rate
|   |-- stress/                  - Robustness tests
|   |   `-- test_stress.py       - GRF completion; high-Mach and low-beta regimes
|   |-- structural/              - Output format tests
|   |   `-- test_structural.py   - HDF5 structure, field names, snapshot counts,
|   |                              thermal pressure identity
|   `-- validation/              - Quantitative validation
|       |-- test_validation.py   - OT vortex vs Stone et al. (2008) published
|       |                          bounds; grid convergence N=64/128/256;
|       |                          checksum guard on reference file
|       |-- athinput.orszag_tang - Athena++ input file for reference generation
|       `-- ot_reference.h5      - Athena++ 500x500 OT vortex reference (t=0.5)
|-- scripts/
|   |-- run_mhd.py               - Run single simulation and animate output
|   |-- run_mhd_MC.py            - Run Monte Carlo ensemble
|   |-- run_tests.py             - Unified test runner (all suites)
|   |-- convert_athena_to_reference.py - Build OT validation reference from
|   |                              an Athena++ VTK output file
|   `-- diagnose_validation.py   - Visual diagnostic for OT validation tests;
|                                  plots CHIMERA vs Athena++ side-by-side
|-- Makefile
`-- README.md
```

---

## Getting Started

### Prerequisites

- Fortran compiler: `gfortran >= 9` or Intel `ifx`/`ifort`
- HDF5 library with Fortran bindings (e.g. `libhdf5-fortran-dev`)
- OpenMP
- Python 3 with `h5py`, `numpy`, `matplotlib`, `pyvista`

### Build

```bash
git clone https://github.com/tjnb0/CHIMERA.git
cd CHIMERA
make          # builds chimera executable
make tests    # also builds Fortran unit test binaries in tests/fortran/
```

Edit the `Makefile` to point to your HDF5 installation if needed.

### Run

```bash
# Single run - Orszag-Tang vortex
./chimera 1 $grid_size 0 ./outputs/ orszag_tang.h5

# MHD Rotor with outflow BCs
./chimera 4 $grid_size 0 ./outputs/ rotor.h5

# Monte Carlo ensemble (100 runs)
for seed in $(seq 1 100); do
    ./chimera 5 $grid_size $seed ./outputs/ mc_run_${seed}.h5
done

# GRF stress-test overrides (args 6 and 7; negative = random sampling)
./chimera 5 $grid_size $seed ./outputs/ out.h5 5.0 -1.0   # target M_s = 5
./chimera 5 $grid_size $seed ./outputs/ out.h5 -1.0 0.1   # target beta = 0.1
```

Or using the Python scripts from the `scripts/` directory:

```bash
python scripts/run_mhd.py          # single run with interactive animation
python scripts/run_mhd_MC.py       # full Monte Carlo ensemble
```

**Command-line arguments:** `problem_type  grid_size  seed  output_path  filename.h5  [target_M_s  target_beta]`

End time, limiter selection, and other scheme options are compile-time parameters
in `mhd_config.f90`.

### Test

```bash
# Default suite: structural + physics + stress (no simulation longer than ~30 s)
python scripts/run_tests.py

# Individual suites
python scripts/run_tests.py --structural   # HDF5 format and Fortran unit tests
python scripts/run_tests.py --physics      # mass conservation, positivity, convergence
python scripts/run_tests.py --stress       # high-Mach and low-beta GRF runs
python scripts/run_tests.py --validation   # OT vortex published bounds and grid convergence
python scripts/run_tests.py --all          # everything
```

The validation suite includes a reference comparison test that requires
`tests/validation/ot_reference.h5`. Generate it from an Athena++ VTK output:

```bash
python scripts/convert_athena_to_reference.py OrszagTang.block0.out1.00001.vtk
```

Use `scripts/diagnose_validation.py` to plot a visual comparison of CHIMERA
versus the Athena++ reference before running the full validation suite.

### Slope Limiter Selection

The slope limiter is selected at build time in `mhd_config.f90`:

```fortran
integer, parameter :: slope_limiter = LIMITER_MC        ! default
! integer, parameter :: slope_limiter = LIMITER_VAN_LEER  ! alternative
```

Change the constant and run `make` to recompile. No other files need editing.

---

## Boundary Conditions

Each side of the domain is assigned independently in the problem setup routine:

```fortran
BC_xlo = BC_PERIODIC   ! left
BC_xhi = BC_OUTFLOW    ! right
BC_ylo = BC_FIXED      ! bottom
BC_yhi = BC_INFLOW     ! top
```

Available types defined in `mhd_config.f90`:

| Constant | Value | Behaviour |
|----------|-------|-----------|
| `BC_PERIODIC` | 1 | Circular wrap |
| `BC_OUTFLOW` | 2 | Zero-gradient outflow |
| `BC_FIXED` | 3 | Prescribed ambient state via `f_ambient` argument |
| `BC_INFLOW` | 4 | Driven inflow (e.g., Coronal Mass Ejections) |

---

## Output Format

Results are written to a single HDF5 file per run:

```
/rho/SNAPSHOT1 ... SNAPSHOTn     - density               [N x N, float64]
/P/SNAPSHOT1   ... SNAPSHOTn     - thermal pressure p    [N x N, float64]
/Bx/SNAPSHOT1  ... SNAPSHOTn     - x magnetic field
/By/SNAPSHOT1  ... SNAPSHOTn     - y magnetic field
/Vx/SNAPSHOT1  ... SNAPSHOTn     - x velocity
/Vy/SNAPSHOT1  ... SNAPSHOTn     - y velocity
/time/sim_time                    - snapshot times        [n, float64]
/parameters/gamma                 - adiabatic index (realized)
/parameters/M_s                   - sonic Mach number (realized)
/parameters/beta                  - plasma beta = 2p/B^2 (realized)
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
- van Leer, B. (1977). *Towards the ultimate conservative difference scheme IV.*
  J. Comput. Phys. 23, 263-275. (monotonized central limiter)
- van Leer, B. (1979). *Towards the ultimate conservative difference scheme.*
  J. Comput. Phys.
- Clain, S., Diot, S. & Loubere, R. (2011). *A high-order finite volume method
  for hyperbolic systems: Multi-dimensional Optimal Order Detection (MOOD).*
  J. Comput. Phys.
- Gardiner, T. A. & Stone, J. M. (2005). *An unsplit Godunov method for ideal
  MHD via constrained transport.* J. Comput. Phys. 205, 509-539.
- Stone, J. M. et al. (2008). *Athena: A new code for astrophysical MHD.*
  ApJS 178, 137-177.

---

## License

MIT License - see [LICENSE](LICENSE) for details.
