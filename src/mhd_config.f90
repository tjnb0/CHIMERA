module mhd_config

    !--------------------------------------------------------------------
    ! Purpose: Global configuration parameters for the 2D MHD simulation
    !          including grid size, scheme options, and constants 
    !--------------------------------------------------------------------

    implicit none

    ! Main Parameters
    integer :: problem_type            ! 1 = OT Vortex; 2 = KH Instab; 3 = Advection; 4 = Rotor; 5 = MC
    integer :: seed                    ! Random seed for Monte Carlo
    character(len=255) :: out_path     ! Output directory
    character(len=255) :: h5_filename  ! Output filename
    integer, parameter :: N = 64       ! Number of points on square grid 
    real(8), parameter :: tEnd = 1.0d0 ! Simulation end time 

    ! Params for HD plots with good time evolution
    ! - OT Vortex: N = 3000, tend = 0.5
    ! - KH Instab: " = 1024,   "  = 1.5
    ! - Advection: " = 1024,   "  = 1.0
    ! - MHD Rotor: " = 2000,   "  = 0.9

    ! Advanced parameters
    logical, parameter :: useSlopeLimiting = .true.  ! Enable/disable slope limiting 
    logical, parameter :: upgrade_2_MOOD = .true.   ! Better reconstruction
    integer, parameter :: BC_PERIODIC = 1
    integer, parameter :: BC_FIXED    = 2
    integer, parameter :: BC_OUTFLOW  = 3
    integer :: BC_x                                  ! Default BC for x
    integer :: BC_y                                  ! Default BC for y

    ! Do not change
    real(8), parameter :: tOut = 0.01d0                 ! Interval for writing data to files
    real(8), parameter :: outInterval = 0.01d0*tEnd     ! Interval to output time stamps
    real(8), parameter :: pi    = 3.1415926535897932d0  ! pi 
    real(8), parameter :: twoPi = 6.2831853071795865d0  ! 2 * pi
    real(8), parameter :: fourPi = 12.566370614359173d0 ! 4 * pi
    real(8), parameter :: boxsize = 1.0d0               ! size of the simulation box 
    real(8), parameter :: courant_fac = 0.3d0           ! Courant–Friedrichs–Lewy safety factor 
    real(8), parameter :: P_floor = 1e-12
    real(8), parameter :: rho_floor = 1e-12

end module mhd_config
