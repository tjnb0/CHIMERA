# CHIMERA
### Constrained-transport, High-resolution Ideal MHD with Ensemble Randomization Approach

A finite-volume solver for the compressible ideal magnetohydrodynamics (MHD) equations in two spatial dimensions, written in modern Fortran. Developed as an independent project.

---

## Overview

CHIMERA evolves the compressible ideal MHD equations in conservation form on a uniform 2D Cartesian grid. The numerical scheme combines second-order MUSCL-Hancock reconstruction with a Rusanov (local Lax-Friedrichs) Riemann solver, and uses constrained transport to preserve ∇·**B** = 0 to machine precision throughout the simulation. OpenMP threading accelerates the reconstruction and slope-limiting passes.

Output is written to HDF5, with each field stored as a sequence of snapshots alongside realised physics diagnostics (γ, Mach number, plasma β).

---

## Governing Equations

CHIMERA advances the compressible MHD system in conservation form:

```
∂ρ/∂t  + ∇·(ρv)                       = 0          (mass)
∂(ρv)/∂t + ∇·(ρvv + P*I − BB)         = 0          (momentum)
∂E/∂t  + ∇·((E + P*)v − B(v·B))       = 0          (energy)
∂B/∂t  − ∇×(v×B)                      = 0          (induction)
```

where `P* = p + |B|²/2` is the total (thermal + magnetic) pressure and `E` is the total energy density (internal + kinetic + magnetic). The gas is closed by an ideal equation of state with adiabatic index γ.

---

## Numerical Methods

| Component | Method |
|-----------|--------|
| Spatial discretisation | Cell-centred finite volume on a uniform Cartesian grid |
| Time integration | Predictor–corrector (MUSCL-Hancock); CFL-limited adaptive timestep |
| Reconstruction | 2nd-order MUSCL with optional MOOD fallback to 1st-order at troubled cells |
| Slope limiting | Van Leer harmonic mean limiter (toggleable) |
| Riemann solver | Local Lax-Friedrichs / Rusanov |
| Divergence control | Constrained transport (CT) on staggered face-centred **B**; ∇·**B** monitored every step |
| Parallelism | OpenMP on reconstruction and slope-limiting loops |

### Time-stepping detail

At each step CHIMERA:
1. Reconstructs primitive variables half a timestep forward (prediction), including advection, pressure-gradient, and Lorentz-force terms.
2. Performs MUSCL spatial reconstruction to obtain left/right states at each cell face (with MOOD check to clamp unphysical overshoots).
3. Evaluates Rusanov fluxes for mass, momentum, energy, and transverse **B**.
4. Updates conserved variables and advances face-centred **B** via constrained transport.

---

## Features

- **MUSCL-Hancock predictor–corrector** — second-order accurate in space and time
- **MOOD reconstruction** — cell-by-cell fallback to first order when reconstructed values exceed the local stencil bounds, preventing spurious oscillations near shocks
- **Van Leer slope limiter** — toggled via `useSlopeLimiting` in `mhd_config.f90`
- **Constrained transport** — staggered face-centred **B** updated via the discrete curl of E_z; div B monitored and printed each timestep
- **Fast magnetosonic CFL condition** — timestep limited by `c_f + |v|` with a Courant factor of 0.4
- **Gaussian random field (GRF) initial conditions** — divergence-free velocity (from stream function) and magnetic field (from vector potential) generated from RBF-smoothed random fields; scaled to target V₀ and B₀ with a fast-speed cap; realised γ, sonic Mach number, and plasma β stored in output
- **HDF5 output** — all primitive fields (ρ, P, Bx, By, Vx, Vy) written as snapshot groups with a shared time array and physics parameters (γ, M_s, β)
- **Multiple boundary condition types** — periodic, fixed, and damped outflow ghost cells implemented
- **OpenMP threading** — slope-limiting and MOOD reconstruction loops parallelised with `!$omp parallel do`

---

## Test Problems

CHIMERA ships with five built-in initial conditions (selected via command-line argument):

| `problem_type` | Problem | Physics tested |
|:-:|---------|----------------|
| 1 | **Orszag–Tang vortex** | 2D MHD turbulence, shock–current-sheet interaction |
| 2 | **Kelvin–Helmholtz instability** | Shear-driven magnetic instability and interface roll-up |
| 3 | **Field-loop advection** | Accuracy of CT scheme; passive advection of a magnetic flux tube |
| 4 | **MHD rotor** | High-density rotating disk; strong discontinuities and torsional Alfvén waves |
| 5 | **Monte Carlo / GRF ensemble** | Statistical studies of compressible MHD turbulence with randomised initial conditions |

Recommended resolutions: Orszag–Tang N=3000, KH N=1024, Advection N=1024, Rotor N=2000.

---

## Code Structure

```
.
├── main.f90                 – Time loop, CFL timestep, I/O scheduling
├── mhd_config.f90           – Global parameters (N, tEnd, CFL, floors, BC flags)
├── mhd_init.f90             – Array allocation, grid setup, all initial conditions
│                              including GRF field generation
├── mhd_change_states.f90    – Primitive ↔ conserved variable conversion
├── mhd_derivatives.f90      – Central-difference gradients, Van Leer slope limiter,
│                              MUSCL/MOOD face reconstruction
├── mhd_flux.f90             – Rusanov flux evaluation, conserved-variable update,
│                              constrained transport
├── mhd_field_ops.f90        – Discrete curl (B from Az), div B diagnostic,
│                              face-to-cell B averaging
├── mhd_ghost_bcs.f90        – Damped outflow ghost-cell filling
└── mhd_write_h5.f90         – HDF5 output (snapshots + physics parameters)
```

---

## Getting Started

### Prerequisites

- Fortran compiler — `gfortran ≥ 9` or Intel `ifx`/`ifort`
- HDF5 library with Fortran bindings (e.g. `libhdf5-fortran-dev` on Debian/Ubuntu)
- OpenMP (included with most compilers; disable with `-fno-openmp` if not needed)

### Build

```bash
git clone https://github.com/tjnb0/CHIMERA.git
cd CHIMERA
make
```

Edit the `Makefile` to point to your HDF5 installation if needed.

### Run

```bash
# Single run — Orszag-Tang vortex, seed 0
./chimera 1 0 ./outputs/ orszag_tang.h5

# Monte Carlo ensemble run (problem type 5, varying seed)
for seed in $(seq 1 100); do
    ./chimera 5 $seed ./outputs/ mc_run_${seed}.h5
done
```

**Command-line arguments:** `problem_type  seed  output_path  filename.h5`

Grid resolution, end time, and scheme options are set as compile-time parameters in `mhd_config.f90`.

---

## Output Format

Results are written to a single HDF5 file per run. The file layout is:

```
/rho/SNAPSHOT1 … SNAPSHOTn     – density         [N × N, float64]
/P/SNAPSHOT1   … SNAPSHOTn     – gas pressure    [N × N, float64]
/Bx/SNAPSHOT1  … SNAPSHOTn     – x magnetic field
/By/SNAPSHOT1  … SNAPSHOTn     – y magnetic field
/Vx/SNAPSHOT1  … SNAPSHOTn     – x velocity
/Vy/SNAPSHOT1  … SNAPSHOTn     – y velocity
/time/sim_time                  – snapshot times  [n, float64]
/parameters/gamma               – adiabatic index (realised)
/parameters/M_s                 – sonic Mach number (realised)
/parameters/beta                – plasma β = 2p/B² (realised)
```

Snapshots are written every `tOut = 0.01` time units. Files can be read with h5py in Python:

```python
import h5py, numpy as np
with h5py.File("orszag_tang.h5", "r") as f:
    rho  = np.array(f["rho/SNAPSHOT10"])
    time = np.array(f["time/sim_time"])
```

---

## References

- Tóth, G. (2000). *The ∇·B = 0 constraint in shock-capturing MHD codes.* J. Comput. Phys.
- Evans, C. R. & Hawley, J. F. (1988). *Simulation of magnetohydrodynamic flows: A constrained transport method.* ApJ.
- van Leer, B. (1979). *Towards the ultimate conservative difference scheme.* J. Comput. Phys.
- Gardiner, T. A. & Stone, J. M. (2005). *An unsplit Godunov method for ideal MHD via constrained transport.* J. Comput. Phys.
- Clain, S., Diot, S. & Loubère, R. (2011). *A high-order finite volume method for hyperbolic systems: Multi-dimensional Optimal Order Detection (MOOD).* J. Comput. Phys.

---

## License

MIT License — see [LICENSE](LICENSE) for details.
