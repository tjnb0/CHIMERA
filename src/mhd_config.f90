module mhd_config

    !--------------------------------------------------------------------
    ! Purpose: Global configuration parameters for the 2D MHD simulation
    !          including grid size, scheme options, and constants 
    !--------------------------------------------------------------------

    implicit none

    ! Main Parameters
    integer :: problem_type            ! 1 = OT Vortex; 2 = KH Instab; 3 = Advection; 4 = Rotor; 5 = MC
    integer :: seed                    ! Random seed for Monte Carlo
    integer :: N                       ! Number of points on square grid
    character(len=255) :: out_path     ! Output directory
    character(len=255) :: h5_filename  ! Output filename
    real(8), parameter :: tEnd = 1.0d0 ! Simulation end time

    ! Params for HD plots with good time evolution
    ! - OT Vortex: N = 3000, tEnd = 0.5
    ! - KH Instab: N = 1024, tEnd = 1.5
    ! - Advection: N = 1024, tEnd = 1.0
    ! - MHD Rotor: N = 2000, tEnd = 0.9


    !--- Advanced parameters ---
    ! Slope limiter selection (no effect when useSlopeLimiting = .false.)
    logical, parameter :: useSlopeLimiting = .true.  ! Enable/disable slope limiting
    integer, parameter :: LIMITER_VAN_LEER = 1                ! Van Leer harmonic mean; robust / more diffusive
    integer, parameter :: LIMITER_MC       = 2                ! Monotonized central; less diffuse
    integer, parameter :: slope_limiter    = LIMITER_MC       ! slope limiter selection

    ! Riemann solver selection
    integer, parameter :: RIEMANN_RUSANOV = 1                 ! local Lax-Friedrichs; robust / more diffusive
    integer, parameter :: RIEMANN_HLLE    = 2                 ! Harten-Lax-van Leer-Einfeldt; less diffusive
    integer, parameter :: RIEMANN_HLLD    = 3                 ! Harten-Lax-van Leer-Discontinuities; resolves MHD waves
    integer, parameter :: riemann_solver  = RIEMANN_RUSANOV   ! Riemann solver selection
    logical, parameter :: fallback_2_MOOD = .true.            ! MOOD reconstruction fallback

    ! Boundary condition type constants
    integer, parameter :: BC_PERIODIC = 1  ! Periodic (wrap-around)
    integer, parameter :: BC_FIXED    = 2  ! Fixed (prescribed) state
    integer, parameter :: BC_OUTFLOW  = 3  ! Zero-gradient outflow
    integer, parameter :: BC_INFLOW   = 4  ! Driven inflow (e.g. CME source)

    ! Per-side BC flags -- set in each problem's setup routine.
    ! xlo = left  (i=1),  xhi = right (i=N)
    ! ylo = bottom(j=1),  yhi = top   (j=N)
    integer :: BC_xlo, BC_xhi
    integer :: BC_ylo, BC_yhi

    ! Do not change
    real(8), parameter :: tOut        = 0.01d0                ! Interval for writing data to files
    real(8), parameter :: outInterval = 0.01d0 * tEnd         ! Interval to output time stamps
    real(8), parameter :: pi    = 3.1415926535897932d0        ! pi
    real(8), parameter :: twoPi = 6.2831853071795865d0        ! 2 * pi
    real(8), parameter :: fourPi = 12.566370614359173d0       ! 4 * pi
    real(8), parameter :: boxsize    = 1.0d0                  ! Size of the simulation box
    real(8), parameter :: courant_fac = 0.3d0                 ! CFL safety factor
    real(8), parameter :: P_floor   = 1.0d-12                 ! Minimum absolute pressure floor
    real(8), parameter :: rho_floor = 1.0d-12                 ! Minimum density floor
    real(8), parameter :: e_floor_frac = 1.0d-10              ! Proportional P floor: p >= e_floor_frac*(KE + MagE)
    real(8), parameter :: cf_max = 1.0d2                      ! Maximumum fast magnetosonic speed
    real(8), parameter :: rho_floor_frac = 1.0d-3             ! proportional density floor
    real(8), parameter :: p_th_mood_frac = 1.0d-4             ! relative thermal-pressure threshold
    real(8), parameter :: p_th_ic_frac = 1.0d-2               ! IC thermal-pressure positivity

    ! Target parameters for GRF 
    ! Negative = sample randomly within bounds
    real(8) :: target_M_s    = -1.0d0   ! Sonic Mach number override in command line
    real(8) :: target_beta   = -1.0d0   ! Plasma beta override in command line
    real(8) :: target_gamma  = -1.0d0   ! Adiabatic index override in command line

end module mhd_config
