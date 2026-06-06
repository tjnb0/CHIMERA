module mhd_init
    use mhd_config
    use mhd_field_ops

    !----------------------------------------------------------------------------
    ! Purpose: Declare and allocate all variables (mesh, time, fluid, field, etc)
    !          One subroutine allocate all the arrays and set up the square grid
    !----------------------------------------------------------------------------

    implicit none

    ! Mesh and time variables (mesh, coords, time steps, etc.)
    integer :: outputCount                           ! Counter for the # of outputs written
    integer :: timeStampCount                        ! Counter for time stamp outputs
    integer :: i, j, io                              ! Loop indices and I/O number
    real(8) :: dx                                    ! Grid cell size (unif. in x and y)
    real(8) :: vol                                   ! Cell volume (dx * dx for 2D)
    real(8) :: t                                     ! Current simulation time
    real(8) :: dt                                    ! Current time step size
    real(8) :: gamma                                 ! Ratio of specific heats for ideal MHD
    real(8), allocatable :: xlin(:)                  ! Cell-centered x-coordinates 
    real(8), allocatable :: xlin_node(:)             ! Node-centered x-coordinates
    real(8), allocatable :: X(:,:), Y(:,:)           ! Cell-centered grid coordinates

    ! Final output storage
    integer :: Nsnap_max = int(tEnd / tOut) + 2
    integer :: Nsnap
    real(8), allocatable :: rho_all(:,:,:), P_all(:,:,:)
    real(8), allocatable :: Bx_all(:,:,:), By_all(:,:,:)
    real(8), allocatable :: Vx_all(:,:,:), Vy_all(:,:,:)
    real(8), allocatable :: time_all(:)
    real(8) :: gamma_actual, beta_actual, M_s_actual    

    ! Fluid and field variables (fluid, Bfield, potential, conserved vars, etc)
    real(8), allocatable :: rho(:,:)                 ! Cell-centered mass density 
    real(8), allocatable :: inv_rho(:,:)             ! 1 div. by rho
    real(8), allocatable :: v_x(:,:)                 ! Face-centered x-velocity 
    real(8), allocatable :: v_y(:,:)                 ! Face-centered y-velocity 
    real(8), allocatable :: vx(:,:)                  ! Cell-centered x-velocity 
    real(8), allocatable :: vy(:,:)                  ! Cell-centered y-velocity 
    real(8), allocatable :: P(:,:)                   ! Cell-centered gas pressure 
    real(8), allocatable :: Bmag(:,:)                ! Cell-centered B-field magnitude
    real(8), allocatable :: Az(:,:)                  ! Magnetic vector potential (z-comp.)
    real(8), allocatable :: b_x(:,:), b_y(:,:)       ! Face-centered B-field components
    real(8), allocatable :: Bx(:,:), By(:,:)         ! Cell-centered B-field components
    real(8), allocatable :: Mass(:,:), Energy(:,:)   ! Conserved vars: mass and energy
    real(8), allocatable :: Momx(:,:), Momy(:,:)     ! Conserved var : momenta

    ! Gradient components (x and y) of primitive vars 
    real(8), allocatable :: rho_dx(:,:), rho_dy(:,:) ! density 
    real(8), allocatable :: vx_dx(:,:), vx_dy(:,:)   ! vx
    real(8), allocatable :: vy_dx(:,:), vy_dy(:,:)   ! vy
    real(8), allocatable :: P_dx(:,:), P_dy(:,:)     ! pressure
    real(8), allocatable :: Bx_dx(:,:), Bx_dy(:,:)   ! Bx
    real(8), allocatable :: By_dx(:,:), By_dy(:,:)   ! By

    ! Limited slopes of primitive vars
    real(8), allocatable :: rho_prime(:,:)           ! density
    real(8), allocatable :: vx_prime(:,:)            ! vx
    real(8), allocatable :: vy_prime(:,:)            ! vy
    real(8), allocatable :: P_prime(:,:)             ! pressure  
    real(8), allocatable :: Bx_prime(:,:)            ! Bx
    real(8), allocatable :: By_prime(:,:)            ! By

    ! Reconstructed primitive boundary values (left/right in x and y)
    real(8), allocatable :: rho_XL(:,:), rho_XR(:,:) ! density (x)
    real(8), allocatable :: rho_YL(:,:), rho_YR(:,:) ! density (y)
    real(8), allocatable :: vx_XL(:,:), vx_XR(:,:)   ! vx (x)
    real(8), allocatable :: vx_YL(:,:), vx_YR(:,:)   ! vy (x)
    real(8), allocatable :: vy_XL(:,:), vy_XR(:,:)   ! vx (y)
    real(8), allocatable :: vy_YL(:,:), vy_YR(:,:)   ! vy (y)
    real(8), allocatable :: P_XL(:,:), P_XR(:,:)     ! pressure (x)
    real(8), allocatable :: P_YL(:,:), P_YR(:,:)     ! pressure (y)
    real(8), allocatable :: Bx_XL(:,:), Bx_XR(:,:)   ! Bx (x)
    real(8), allocatable :: Bx_YL(:,:), Bx_YR(:,:)   ! By (x)
    real(8), allocatable :: By_XL(:,:), By_XR(:,:)   ! Bx (y)
    real(8), allocatable :: By_YL(:,:), By_YR(:,:)   ! By (y)

    ! Flux arrays for conserved vars (x and y)
    real(8), allocatable :: flux_Mass_X(:,:)         ! mass through X face
    real(8), allocatable :: flux_Mass_Y(:,:)         ! mass through Y face
    real(8), allocatable :: flux_Momx_X(:,:)         ! x momentum through X face
    real(8), allocatable :: flux_Momx_Y(:,:)         ! x momentum through Y face
    real(8), allocatable :: flux_Momy_X(:,:)         ! y momentum through X face
    real(8), allocatable :: flux_Momy_Y(:,:)         ! y momentum through Y face
    real(8), allocatable :: flux_Energy_X(:,:)       ! energy through X face
    real(8), allocatable :: flux_Energy_Y(:,:)       ! energy through Y face
    real(8), allocatable :: flux_By_X(:,:)           ! By through X face
    real(8), allocatable :: flux_Bx_Y(:,:)           ! Bx through Y face

    ! Wave speeds and divergence
    real(8), allocatable :: c0_sq(:,:)               ! Local sound speed (squared)
    real(8), allocatable :: ca_sq(:,:)               ! Local Alfven speed (squared)
    real(8), allocatable :: cf(:,:)                  ! Local fast magnetosonic speed
    real(8), allocatable :: divB(:,:)                ! Magnetic field divergence

    ! Ghost-cell padded arrays (N+2) x (N+2) for open/outflow boundary conditions.
    ! Interior [2:N+1, 2:N+1] holds the physical domain; rows/columns 1 and N+2
    ! are ghost cells filled each timestep by fill_ghost_cells() in mhd_bc.
    real(8), allocatable :: rho_pad(:,:), vx_pad(:,:), vy_pad(:,:)
    real(8), allocatable :: P_pad(:,:),   Bx_pad(:,:), By_pad(:,:)

    
contains


    subroutine initialize_variables()
        ! 
        ! Allocate all arrays in sim. All 1D are size N, all 2D are size NxN 
        ! 

        ! Allocate 1D coordinate arrays for mesh (cell centers and nodes)
        allocate(xlin(N))

        ! Allocate 2D grid arrays for cell centers
        allocate(X(N,N), Y(N,N))

        ! Allocate magnetic field and potential arrays
        allocate(Az(N,N))               ! Magnetic vector potential (z-component)
        allocate(b_x(N,N), b_y(N,N))    ! Face-centered B-field 
        allocate(Bx(N,N), By(N,N))      ! Cell-centered B-field 

        ! Allocate conserved variable arrays (mass, momentum, energy)
        allocate(Mass(N,N), Momx(N,N), Momy(N,N), Energy(N,N))

        ! Allocate primitive variable arrays (cell- and face-centered)
        allocate(rho(N,N), vx(N,N), vy(N,N), P(N,N), Bmag(N,N))
        allocate(v_x(N,N), v_y(N,N))

        ! Allocate gradient arrays for primitive vars
        allocate(rho_dx(N,N), rho_dy(N,N))    ! Density gradients
        allocate(vx_dx(N,N), vx_dy(N,N))      ! x-velocity gradients
        allocate(vy_dx(N,N), vy_dy(N,N))      ! y-velocity gradients
        allocate(P_dx(N,N), P_dy(N,N))        ! Pressure gradients
        allocate(Bx_dx(N,N), Bx_dy(N,N))      ! Bx gradients
        allocate(By_dx(N,N), By_dy(N,N))      ! By gradients

        ! Allocate arrays for limited slopes (after slope limiting)
        allocate(rho_prime(N,N), vx_prime(N,N), vy_prime(N,N), P_prime(N,N))
        allocate(Bx_prime(N,N), By_prime(N,N))

        ! Allocate arrays for reconstructed boundary values (left/right in x and y)
        allocate(rho_XL(N,N), rho_XR(N,N), rho_YL(N,N), rho_YR(N,N))
        allocate(vx_XL(N,N), vx_XR(N,N), vx_YL(N,N), vx_YR(N,N))
        allocate(vy_XL(N,N), vy_XR(N,N), vy_YL(N,N), vy_YR(N,N))
        allocate(P_XL(N,N), P_XR(N,N), P_YL(N,N), P_YR(N,N))
        allocate(Bx_XL(N,N), Bx_XR(N,N), Bx_YL(N,N), Bx_YR(N,N))
        allocate(By_XL(N,N), By_XR(N,N), By_YL(N,N), By_YR(N,N))

        ! Allocate flux arrays for all conserved variables in both x and y directions
        allocate(flux_Mass_X(N,N), flux_Mass_Y(N,N))
        allocate(flux_Momx_X(N,N), flux_Momx_Y(N,N))
        allocate(flux_Momy_X(N,N), flux_Momy_Y(N,N))
        allocate(flux_Energy_X(N,N), flux_Energy_Y(N,N))
        allocate(flux_By_X(N,N), flux_Bx_Y(N,N))

        ! Allocate arrays for wave speeds and divergence
        allocate(c0_sq(N,N))  ! Local sound speed
        allocate(ca_sq(N,N))  ! Local Alfven speed
        allocate(cf(N,N))     ! Local fast magnetosonic speed
        allocate(divB(N,N))   ! Magnetic field divergence

        ! Allocate ghost-cell padded arrays for open boundary conditions
        allocate(rho_pad(N+2,N+2), vx_pad(N+2,N+2), vy_pad(N+2,N+2))
        allocate(P_pad(N+2,N+2),   Bx_pad(N+2,N+2), By_pad(N+2,N+2))

        ! Allocate final output arrays
        allocate(rho_all(Nsnap_max, N, N))
        allocate(P_all(Nsnap_max, N, N))
        allocate(Bx_all(Nsnap_max, N, N))
        allocate(By_all(Nsnap_max, N, N))
        allocate(Vx_all(Nsnap_max, N, N))
        allocate(Vy_all(Nsnap_max, N, N))
        allocate(time_all(Nsnap_max))

    end subroutine initialize_variables


    subroutine get_user_options()
    !
    !   Parse command-line arguments.
    !
    !   Usage: ./mhd_sim <problem_type> <N> <seed> [output_path] [h5_filename]
    !
    !   Required:
    !       problem_type : 1=Orszag-Tang, 2=KH, 3=Field Loop, 4=Rotor, 5=MC
    !       N            : grid size (N x N cells, must be >= 4)
    !       seed         : random seed for Monte Carlo (ignored for problems 1-4)
    !
    !   Optional:
    !       output_path  : default './outputs/'
    !       h5_filename  : default 'primitive_snaps.h5'
    !
        character(len=256) :: arg_str
        integer :: nargs, istat

        nargs = command_argument_count()

        if (nargs < 3) then
            print *, "ERROR: Usage: ./mhd_sim <problem_type> <N> <seed> [output_path] [h5_filename]"
            print *, "  problem_type : 1=Orszag-Tang 2=KH 3=Field-Loop 4=Rotor 5=MC"
            print *, "  N            : grid size (integer >= 4)"
            print *, "  seed         : random seed (integer)"
            print *, "  output_path  : optional (default: './outputs/')"
            print *, "  h5_filename  : optional (default: 'primitive_snaps.h5')"
            stop
        end if

        ! Arg 1: problem_type
        call get_command_argument(1, arg_str)
        read(arg_str, *, iostat=istat) problem_type
        if (istat /= 0 .or. problem_type < 1 .or. problem_type > 5) then
            print *, "ERROR: problem_type must be 1-5, got: '", trim(arg_str), "'"
            stop
        end if

        ! Arg 2: N (grid size)
        call get_command_argument(2, arg_str)
        read(arg_str, *, iostat=istat) N
        if (istat /= 0 .or. N < 4) then
            print *, "ERROR: N must be an integer >= 4, got: '", trim(arg_str), "'"
            stop
        end if

        ! Arg 3: seed
        call get_command_argument(3, arg_str)
        read(arg_str, *, iostat=istat) seed
        if (istat /= 0) then
            print *, "ERROR: Invalid seed: '", trim(arg_str), "'"
            stop
        end if

        ! Arg 4: output path (optional)
        if (nargs >= 4) then
            call get_command_argument(4, out_path)
            out_path = trim(adjustl(out_path))
        else
            out_path = "./outputs/"
        end if

        ! Arg 5: HDF5 filename (optional)
        if (nargs >= 5) then
            call get_command_argument(5, h5_filename)
        else
            h5_filename = "primitive_snaps.h5"
        end if

    end subroutine get_user_options


    subroutine initialize_problem()
    !
    !   Setup the initial conditions for the desired problem
    !
        select case (problem_type)
            case (1)
                gamma = 5.0d0/3.0d0 ! Specific heats ratio
                call setup_OT_vortex()

            case (2)
                gamma = 1.4d0
                call setup_KH_instab()

            case (3)
                gamma = 5.0d0/3.0d0
                call setup_FL_advec()

            case(4)
                gamma = 1.4d0
                call setup_rotor()
                     
            case(5)
                call setup_Monte_Carlo()
            
            case default
                print *, 'ERROR: Unrecognized problem_type in initialize_problem'
                stop
        end select

    end subroutine initialize_problem


    !--------------------------
    !--- Orszag-Tang vortex ---
    !--------------------------
    ! From: ATHENA - Princeton
    !
    subroutine setup_OT_vortex()   
        
        BC_xlo = BC_PERIODIC;  BC_xhi = BC_PERIODIC
        BC_ylo = BC_PERIODIC;  BC_yhi = BC_PERIODIC

        rho = 25.0d0 / (36.0d0 * pi) ! uniform density
        P   =  5.0d0 / (12.0d0 * pi) ! uniform gas pressure
        vx  = -sin(twopi * Y)        ! inital x velocity 
        vy  =  sin(twopi * X)        ! inital y velocity
        Az = cos(fourPi * X) / (fourPi * sqrt(fourPi)) + &
                cos(twopi  * Y) / (twopi * sqrt(fourPi))

        call compute_curl_2d(Az, dx, N, N, b_x, b_y)
        call average_face_to_cell_B(b_x, b_y, N, N, Bx, By)
        P = P + 0.5d0 * (Bx*Bx + By*By)
    end subroutine setup_OT_vortex


    !------------------------------------
    !--- Kelvin-Helmholtz Instability ---
    !------------------------------------
    ! From: ATHENA - Princeton
    !
    subroutine setup_KH_instab()
        integer :: ix, iy
        real(8) :: x0, y0, sigma, amp
        
        BC_xlo = BC_PERIODIC;  BC_xhi = BC_PERIODIC
        BC_ylo = BC_PERIODIC;  BC_yhi = BC_PERIODIC

        P  = 2.5d0    ! Constant pressure
        Bx = 0.0d0    ! Uniform magnetic field
        By = 0.0d0    ! Uniform magnetic field
        do iy = 1, N
            do ix = 1, N
                ! Density and x-velocity shear
                if (abs(Y(ix, iy) - boxsize/2.d0) <= 0.25d0 * boxsize) then
                    rho(ix, iy) = 2.0d0
                    vx(ix, iy)  = 0.5d0
                else
                    rho(ix, iy) = 1.0d0
                    vx(ix, iy)  = -0.5d0
                end if

                ! Single Gaussian perturbation at middle of lower interface 
                sigma = 0.05d0 * boxsize  ! Width of the perturbation
                x0 = boxsize / 2.d0       ! Center at x = boxsize/2
                y0 = 0.25d0 * boxsize     ! Lower interface
                amp = 0.5d0               ! Amplitude of perturbation 
                vy(ix, iy) = amp * exp(-((X(ix, iy)-x0)**2 + (Y(ix, iy)-y0)**2) &
                            / (2.d0 * sigma**2))
            end do
        end do
        if (allocated(Bmag)) Bmag = sqrt(Bx**2 + By**2)
        if (allocated(Az))   Az = 0.d0
        b_x = 0.0d0;  b_y = 0.0d0       ! face-centred B: B=0 everywhere for KH
    end subroutine setup_KH_instab


    !----------------------------
    !--- Field Loop Advection ---
    !----------------------------
    ! From: FLASH - Rochester
    !
    subroutine setup_FL_advec()
        integer :: ix, iy
        real(8) :: x0, y0, r0, r, A0

        BC_xlo = BC_PERIODIC;  BC_xhi = BC_PERIODIC
        BC_ylo = BC_PERIODIC;  BC_yhi = BC_PERIODIC

        r0 = 0.3d0 * boxsize ! Field loop radius
        A0 = 1.0d-3          ! Field loop amplitude
        rho = 1.0d0;        P  = 1.0d0
        vx  = -2.0d0;        vy = -1.0d0
        x0 = boxsize/2.0d0; y0 = boxsize/2.0d0
        do iy = 1, N
            do ix = 1, N
                ! Build vector potential field Az for centered difference later
                r = sqrt( (X(ix, iy) - x0)**2 + (Y(ix, iy) - y0)**2 )
                if (r < r0) then
                    Az(ix, iy) = A0 * (r0 - r)
                else
                    Az(ix, iy) = 0.0d0
                end if
            end do
        end do

        ! Compute B = curl(Az)
        call compute_curl_2d(Az, dx, N, N, b_x, b_y)
        call average_face_to_cell_B(b_x, b_y, N, N, Bx, By)
        P = P + 0.5d0 * (Bx*Bx + By*By)

        ! Set magnetic field magnitude array
        if (allocated(Bmag)) Bmag = sqrt(Bx**2 + By**2)
    end subroutine setup_FL_advec
    

    !-------------------
    !--- MHD Rotor -----
    !-------------------
    ! From: FLASH - Rochester
    !
    subroutine setup_rotor()
        integer :: ix, iy
        real(8) :: x0, y0, r0, r, rx, ry, r1, u0, f

        BC_xlo = BC_OUTFLOW;  BC_xhi = BC_OUTFLOW
        BC_ylo = BC_OUTFLOW;  BC_yhi = BC_OUTFLOW

        rho = 1.0d0
        x0    = boxsize/2.d0;          y0  = boxsize/2.d0
        r0    = 0.1d0  * boxsize;      r1  = 0.115d0 * boxsize
        u0    = 2.0d0;                 P   = 1.0d0            
        Bx    = 5.0d0 / sqrt(4.d0*pi); By  = 0.0d0 
        vx    = 0.0d0;                 vy  = 0.0d0

        do iy = 1, N
            do ix = 1, N
                rx = X(ix,iy) - x0
                ry = Y(ix,iy) - y0
                r  = sqrt(rx*rx + ry*ry)

                if (r <= r0) then
                    ! Core: constant density, solid-body rotation
                    rho(ix,iy) = 10.0d0
                    if (r0 > 0.d0) then
                        vx(ix,iy) = -u0 * (ry / r0)
                        vy(ix,iy) =  u0 * (rx / r0)
                    else
                        vx(ix,iy) = 0.d0
                        vy(ix,iy) = 0.d0
                    end if

                else if (r < r1) then
                    ! Linear taper between r0 and r1
                    f = (r1 - r) / (r1 - r0)   ! f ->1 at r=r0, ->0 at r=r1
                    rho(ix,iy) = 1.0d0 + 9.0d0 * f
                    ! Taper both magnitude and 1/r to preserve angular profile
                    if (r > 0.d0) then
                        vx(ix,iy) = -f * u0 * (ry / r)
                        vy(ix,iy) =  f * u0 * (rx / r)
                    else
                        vx(ix,iy) = 0.d0
                        vy(ix,iy) = 0.d0
                    end if
                end if
            end do
        end do
        P = P + 0.5d0 * (Bx*Bx + By*By)

        if (allocated(Bmag)) Bmag = sqrt(Bx**2 + By**2)
        if (allocated(Az))   Az   = 0.d0
        b_x = Bx;  b_y = 0.0d0          ! face-centred B must match cell-centred initial state
    end subroutine setup_rotor


    !-------------------
    !--- Monte Carlo ---
    !-------------------
    ! Produces a randomized sample based on parameterizing mach number (M_s),
    ! plasma beta, and gamma. These are uniformly sampled from predefined bounds.
    ! Uses Gaussian Random Field (GRF) approach following Rosofsky & Huerta (2023)
    ! to generate smooth, divergence-free initial conditions.
    !
    subroutine setup_Monte_Carlo()
        integer :: rng_size
        integer, allocatable :: rng_seed(:)
        
        ! Physics parameters (to be sampled and passed to setup routines)
        real(8) :: gamma_val, M_s, beta_val
        real(8) :: V_0, B_0, c_s, V_A, M_A

        ! Initialize RNG with seed
        call random_seed(size=rng_size)
        allocate(rng_seed(rng_size))
        rng_seed(1) = seed
        do i = 2, rng_size
            rng_seed(i) = 1103515245 * rng_seed(i-1) + 12345
            rng_seed(i) = mod(rng_seed(i), 2**30)
        end do
        call random_seed(put=rng_seed)
        deallocate(rng_seed)

        ! Sample physics parameters: (gamma, M_s, beta) -> (V_0, B_0, c_s, ...)
        call sample_physics_parameters(gamma_val, M_s, beta_val, V_0, B_0, c_s, V_A, M_A)

        ! Generate initial conditions via Gaussian Random Fields
        call setup_GRF_fields(gamma_val, V_0, B_0)

    end subroutine setup_Monte_Carlo


    !-------------------------------------
    !--- Sample Physics Parameters -------
    !-------------------------------------
    ! Samples fundamental dimensionless parameters and computes derived scales
    ! 
    ! Outputs:
    !   - gamma_val : Specific heat ratio
    !   - M_s       : Sonic Mach number (V_0 / c_s)
    !   - beta_val  : Plasma beta (2*P_gas / B^2)
    !   - V_0       : Velocity scale [derived]
    !   - B_0       : Magnetic field scale [derived]
    !   - c_s       : Sound speed [derived]
    !   - V_A       : Alfven speed [derived]
    !   - M_A       : Alfven Mach number [derived]
    !
    subroutine sample_physics_parameters(gamma_val, M_s, beta_val, V_0, B_0, c_s, V_A, M_A)
        real(8), intent(out) :: gamma_val, M_s, beta_val
        real(8), intent(out) :: V_0, B_0, c_s, V_A, M_A
        
        real(8) :: rho_mean, P_mean
        real(8) :: rand_vals(3)
        

        ! STEP 1: Sample primary dimensionless parameters
        !
        call random_number(rand_vals)
        
        ! Gamma: specific heat ratio
        gamma_val = 1.2d0 + 0.5d0 * rand_vals(1)  ! Range: Unif(1.2, 1.7)
        !   - gamma ~ 5/3 ~ 1.67: Monatomic gas
        !   - gamma ~ 1.4: Diatomic gas
        !   - gamma ~ 1.2: Approaching isothermal
        
        ! Sonic Mach number: M_s = V_0 / c_s
        M_s = 10.0d0**(log10(0.2d0) + (log10(2.0d0) - log10(0.2d0)) * rand_vals(2))
        M_s = min(max(M_s, 0.2d0), 2.0d0)
        ! Sample log10(M_s) uniformly in [log10(0.1), log10(3.0)]
        !   - ~33% in [0.2, 0.55)  : Subsonic
        !   - ~33% in [0.55, 1.7)  : Transonic/sonic
        !   - ~33% in [1.7, 2.0]   : Supersonic
        
        ! Plasma beta: beta = 2*P_gas / B^2 
        beta_val = 10.0d0**(log10(0.5d0) + (log10(5.0d0) - log10(0.5d0)) * rand_vals(3))
        beta_val = min(max(beta_val, 0.5d0), 5.0d0)
        ! Sample log10(beta) uniformly in [log10(0.1), log10(5.0)]
        !   - ~33% in [0.5, 0.7)   : Moderately to strongly magnetized
        !   - ~33% in [0.7, 2.2)   : Moderate magnetization (near equipartition)
        !   - ~33% in [2.2, 5.0]   : Weakly magnetized (gas pressure dominated)
        
        ! STEP 2: Fix normalized base state
        rho_mean = 1.0d0
        P_mean   = 1.0d0
        
        ! STEP 3: Derive physical scales from dimensionless parameters
        !
        ! Sound speed: c_s = sqrt(gamma * P / rho)
        c_s = sqrt(gamma_val * P_mean / rho_mean)
        
        ! Velocity scale: V_0 = M_s * c_s
        V_0 = M_s * c_s
        
        ! Magnetic field scale: B_0 = sqrt(2 * P / beta)
        B_0 = sqrt(2.0d0 * P_mean / beta_val)
        
        ! Alfven speed: V_A = B_0 / sqrt(rho)
        V_A = B_0 / sqrt(rho_mean)
        
        ! Alfven Mach number (derived): M_A = V_0 / V_A = M_s * sqrt(gamma * beta / 2)
        M_A = V_0 / (V_A + 1.0d-10)  ! Avoid division by zero

        ! STEP 4: Print sampled parameters        
        !print *, "=========================================="
        !print *, "Sampled Physics Parameters"
        !print *, "=========================================="
        !print *, "Primary Parameters (sampled):"
        !print *, "  gamma            =", gamma_val
        !print *, "  M_s (sonic Mach) =", M_s
        !print *, "  beta (plasma)    =", beta_val
        !print *, ""
        !print *, "Derived Scales:"
        !print *, "  c_s (sound speed)    =", c_s
        !print *, "  V_0 (velocity scale) =", V_0
        !print *, "  B_0 (B-field scale)  =", B_0
        !print *, "  V_A (Alfven speed)   =", V_A
        !print *, "  M_A (Alfven Mach)    =", M_A
        !print *, "=========================================="
        
    end subroutine sample_physics_parameters


    !---------------------------------------
    !--- Gaussian Random Field Generator ---
    !---------------------------------------
    ! Generates a single Gaussian Random Field on [0, boxsize]^2 using the
    ! radial basis function (RBF) kernel in Fourier space to enforce periodic BCs.
    !
    ! Method follows Rosofsky & Huerta (2023), Sec. 4.1:
    !   k_l(x1, x2) = exp(-||x1 - x2||^2 / (2 * l^2))
    !
    ! In Fourier space the RBF kernel has a Gaussian power spectrum:
    !   S(k) = exp(-0.5 * l^2 * |k|^2)
    ! where k is the physical wavenumber vector (rad / length unit).
    !
    ! The GRF is sampled by:
    !   1. Drawing independent complex Gaussian noise at each (kx, ky).
    !   2. Multiplying by sqrt(S(k)) to impose the RBF power spectrum.
    !   3. Inverse DFT back to real space (only the real part is kept).
    !
    ! Implementation note:
    !   The IDFT is computed as a sum of cosine/sine modes rather than
    !   calling an external FFT library, keeping the module self-contained.
    !   For grids with N > ~64 consider replacing the inner loops with a
    !   call to an FFT library (e.g. FFTW or MKL) for O(N^2 log N) scaling.
    !
    ! Inputs:
    !   l_scale  : RBF kernel length scale (controls spatial smoothness)
    !   field    : Output NxN real-valued GRF
    !
    subroutine generate_GRF(l_scale, field)
        real(8), intent(in)  :: l_scale
        real(8), intent(out) :: field(N, N)

        ! Fourier-space coefficient arrays (real and imaginary parts)
        real(8), allocatable :: fhat_re(:,:), fhat_im(:,:)

        ! Local scalars - all declared before any executable statements
        real(8) :: kx_val, ky_val, k_sq, S_k
        real(8) :: rand_u1, rand_u2, gauss_re, gauss_im
        real(8) :: phase_arg, norm_re
        integer :: ix, iy, kx_idx, ky_idx, kx_signed, ky_signed

        allocate(fhat_re(N, N), fhat_im(N, N))

        ! ----------------------------------------------------------------
        ! Step 1: Build spectral coefficients fhat(kx, ky)
        !
        ! Array layout: index 1..N maps to signed wavenumber 0..N/2-1, -N/2..-1
        ! (standard DFT ordering). Physical wavenumber:
        !   kx_phys = 2*pi * kx_signed / boxsize
        ! ----------------------------------------------------------------
        do ky_idx = 1, N
            do kx_idx = 1, N

                kx_signed = kx_idx - 1
                ky_signed = ky_idx - 1
                if (kx_signed > N/2) kx_signed = kx_signed - N
                if (ky_signed > N/2) ky_signed = ky_signed - N

                kx_val = twopi * dble(kx_signed) / boxsize
                ky_val = twopi * dble(ky_signed) / boxsize
                k_sq   = kx_val**2 + ky_val**2

                ! sqrt(S(k)): scale amplitudes so power ~ exp(-l^2*k^2/2)
                S_k = exp(-0.5d0 * l_scale**2 * k_sq)

                ! Box-Muller: two uniform samples -> complex Gaussian N(0,1)+iN(0,1)
                call random_number(rand_u1)
                call random_number(rand_u2)
                rand_u1  = max(rand_u1, 1.0d-15)
                gauss_re = sqrt(-2.0d0 * log(rand_u1)) * cos(twopi * rand_u2)
                gauss_im = sqrt(-2.0d0 * log(rand_u1)) * sin(twopi * rand_u2)

                fhat_re(kx_idx, ky_idx) = S_k * gauss_re
                fhat_im(kx_idx, ky_idx) = S_k * gauss_im

            end do
        end do

        ! ----------------------------------------------------------------
        ! Step 2: Enforce Hermitian symmetry so IDFT gives a real field.
        ! Zero and Nyquist modes must be purely real.
        ! ----------------------------------------------------------------
        fhat_im(1,     1)     = 0.0d0
        fhat_im(N/2+1, 1)     = 0.0d0
        fhat_im(1,     N/2+1) = 0.0d0
        fhat_im(N/2+1, N/2+1) = 0.0d0

        ! ----------------------------------------------------------------
        ! Step 3: Inverse DFT -> real-space field
        !
        ! f(ix,iy) = (1/N^2) * Re[ sum_{kx,ky} fhat(kx,ky)
        !                          * exp(i*2*pi*((kx-1)*(ix-1)+(ky-1)*(iy-1))/N) ]
        !
        ! NOTE: This is O(N^4). For large N replace with an FFT call, e.g.:
        !   fhat_cmplx = cmplx(fhat_re, fhat_im)
        !   call ifft2d(fhat_cmplx, field)   ! via FFTW dzfft2d / c_dfftw_plan
        ! ----------------------------------------------------------------
        field = 0.0d0
        do iy = 1, N
            do ix = 1, N
                norm_re = 0.0d0
                do ky_idx = 1, N
                    do kx_idx = 1, N
                        phase_arg = twopi * ( dble((kx_idx-1)*(ix-1)) &
                                            + dble((ky_idx-1)*(iy-1)) ) / dble(N)
                        norm_re = norm_re &
                                + fhat_re(kx_idx, ky_idx) * cos(phase_arg) &
                                - fhat_im(kx_idx, ky_idx) * sin(phase_arg)
                    end do
                end do
                field(ix, iy) = norm_re / dble(N*N)
            end do
        end do

        ! ----------------------------------------------------------------
        ! Step 4: Remove DC offset so the field has zero mean
        ! ----------------------------------------------------------------
        field = field - sum(field) / dble(N*N)

        deallocate(fhat_re, fhat_im)

    end subroutine generate_GRF


    !----------------------------------------------
    !--- Gaussian Random Field Initial Conditions --
    !----------------------------------------------
    ! Generates smooth, divergence-free MHD initial conditions following
    ! Rosofsky & Huerta (2023), Sec. 4.1.
    !
    ! Procedure:
    !   1. Sample two independent GRFs using the RBF kernel with l = l_scale
    !   2. Use GRF_1 as the vorticity stream function psi:
    !        v = c_psi * curl(psi) = c_psi * (d_psi/dy, -d_psi/dx)
    !      This guarantees div(v) = 0 identically.
    !   3. Use GRF_2 as the magnetic vector potential Az:
    !        B = c_A * curl(Az) = c_A * (dAz/dy, -dAz/dx)
    !      This guarantees div(B) = 0 identically.
    !   4. Scale v and B by c_psi = 0.1 and c_A = 0.005 (Rosofsky & Huerta values),
    !      then additionally rescale to match target V_0 and B_0 from physics sampling.
    !   5. Set uniform thermodynamic base state with small GRF perturbations.
    !
    subroutine setup_GRF_fields(gamma_val, V_0, B_0)
        real(8), intent(in) :: gamma_val, V_0, B_0

        real(8), allocatable :: psi(:,:), A_pot(:,:), grf_rho(:,:), grf_P(:,:)

        real(8), parameter :: c_psi = 0.05d0
        real(8), parameter :: c_A   = 0.001d0
        real(8), parameter :: rho_mean      = 1.0d0
        real(8), parameter :: P_mean        = 1.0d0
        real(8), parameter :: rho_pert_amp  = 0.1d0
        real(8), parameter :: P_pert_amp    = 0.05d0
        real(8), parameter :: target_cf_max = 3.0d0   ! max fast speed (tune with courant_fac)

        real(8) :: v_rms, B_rms, c_s_avg, P_avg
        real(8) :: c_s_loc, c_a_loc, c_f_loc, v_loc, total_max, scale
        real(8) :: l_scale, runif
        integer :: ix, iy, ixp, ixm, iyp, iym

        BC_xlo = BC_PERIODIC;  BC_xhi = BC_PERIODIC
        BC_ylo = BC_PERIODIC;  BC_yhi = BC_PERIODIC
        gamma = gamma_val

        allocate(psi(N,N), A_pot(N,N), grf_rho(N,N), grf_P(N,N))

        ! Step 1: Generate GRFs
        call random_number(runif)
        l_scale = 0.15d0 + 0.1d0 * runif
        call generate_GRF(l_scale, psi)
        call generate_GRF(l_scale, A_pot)
        call generate_GRF(l_scale, grf_rho)
        call generate_GRF(l_scale, grf_P)

        ! Step 2: Divergence-free velocity from stream function psi
        do iy = 1, N
            do ix = 1, N
                ixp = mod(ix,    N) + 1
                ixm = mod(ix-2+N, N) + 1
                iyp = mod(iy,    N) + 1
                iym = mod(iy-2+N, N) + 1
                vx(ix,iy) =  c_psi * (psi(ix,iyp) - psi(ix,iym)) / (2.0d0*dx)
                vy(ix,iy) = -c_psi * (psi(ixp,iy) - psi(ixm,iy)) / (2.0d0*dx)
            end do
        end do

        ! Step 3: Scale velocity to target V_0
        v_rms = sqrt(sum(vx**2 + vy**2) / dble(N*N))
        if (v_rms > 1.0d-10) then
            vx = vx * (V_0 / v_rms)
            vy = vy * (V_0 / v_rms)
        end if

        ! Step 4: Divergence-free B from magnetic potential
        Az = c_A * A_pot
        call compute_curl_2d(Az, dx, N, N, b_x, b_y)
        call average_face_to_cell_B(b_x, b_y, N, N, Bx, By)

        ! Step 5: Scale B to target B_0
        B_rms = sqrt(sum(Bx**2 + By**2) / dble(N*N))
        if (B_rms > 1.0d-10) then
            Az = Az * (B_0 / B_rms)
            Bx = Bx * (B_0 / B_rms)
            By = By * (B_0 / B_rms)
            ! Also rescale face-centered fields so CT stays consistent
            b_x = b_x * (B_0 / B_rms)
            b_y = b_y * (B_0 / B_rms)
        else
            Az = 0.0d0
            Bx = 0.1d0 * B_0
            By = 0.0d0
        end if
        if (allocated(Bmag)) Bmag = sqrt(Bx**2 + By**2)

        ! Step 6: Set rho and P BEFORE the speed cap (they're needed for c_s)
        if (maxval(abs(grf_rho)) > 1.0d-10) &
            grf_rho = grf_rho / maxval(abs(grf_rho))
        if (maxval(abs(grf_P)) > 1.0d-10) &
            grf_P = grf_P / maxval(abs(grf_P))

        rho = max(rho_mean * (1.0d0 + rho_pert_amp * grf_rho), 0.1d0)
        P   = max(P_mean   * (1.0d0 + P_pert_amp   * grf_P),   0.1d0)

        ! Step 7: Speed cap - only rescale VELOCITY, not B
        !   B contributes to c_f but is already set to the physically correct
        !   scale via beta. Squashing it here would break the sampled beta.
        total_max = 0.0d0
        do iy = 1, N
            do ix = 1, N
                c_s_loc   = sqrt(gamma * P(ix,iy) / rho(ix,iy))
                c_a_loc   = sqrt((Bx(ix,iy)**2 + By(ix,iy)**2) / rho(ix,iy))
                c_f_loc   = sqrt(c_s_loc**2 + c_a_loc**2)
                v_loc     = sqrt(vx(ix,iy)**2 + vy(ix,iy)**2)
                total_max = max(total_max, c_f_loc + v_loc)
            end do
        end do

        if (total_max > target_cf_max) then
            scale = target_cf_max / total_max
            vx = vx * scale
            vy = vy * scale
            ! B is NOT rescaled here - beta is a sampled physics parameter
        end if

        ! Step 8: Diagnostics
        v_rms   = sqrt(sum(vx**2 + vy**2) / dble(N*N))
        B_rms   = sqrt(sum(Bx**2 + By**2) / dble(N*N))
        P_avg   = sum(P) / dble(N*N)
        c_s_avg = sqrt(gamma * P_avg / (sum(rho) / dble(N*N)))

        gamma_actual = gamma
        M_s_actual   = v_rms  / (c_s_avg + 1.0d-16)
        beta_actual  = 2.0d0 * P_avg / (B_rms**2 + 1.0d-16)

        !print *, ""
        !print *, "GRF Initial Conditions - Realized Values:"
        !print *, "  gamma =", gamma_actual
        !print *, "  M_s   =", M_s_actual
        !print *, "  beta  =", beta_actual
        !print *, "  l_scale (RBF) =", l_scale
        !print *, ""
        !print *, "Field Statistics:"
        !print *, "  rho: min/mean/max =", minval(rho), sum(rho)/dble(N*N), maxval(rho)
        !print *, "  |v|: RMS/max      =", v_rms, maxval(sqrt(vx**2 + vy**2))
        !print *, "  |B|: RMS/max      =", B_rms, maxval(Bmag)
        !print *, ""

        ! Convert to total pressure for setup
        P = P + 0.5d0*(Bx*Bx + By*By)
        deallocate(psi, A_pot, grf_rho, grf_P)

    end subroutine setup_GRF_fields


end module mhd_init