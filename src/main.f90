program main
    use mhd_config
    use mhd_init
    use mhd_field_ops
    use mhd_change_states
    use mhd_derivatives
    use mhd_flux
    use mhd_write_h5
    use mhd_bc

    implicit none
    ! For trimmed tensors
    integer :: actual_nsnap
    real(kind=8), dimension(:,:,:), allocatable :: rho_trim, P_trim, Bx_trim, By_trim, Vx_trim, Vy_trim
    real(kind=8), dimension(:), allocatable :: time_trim
    character(len=512) :: full_h5_path


    ! Get desired problem, output path, and filename
    call get_user_options()

    ! Setup simulation
    call initialize_variables()             ! Initialize arrays and variables

    ! Generate uniform grid
    dx   = boxsize / dble(N)                ! cell width
    vol  = dx*dx                            ! cell area
    xlin = [( (i - 0.5d0), i = 1, N)] * dx  ! cell centers

    ! Fill 2D grid cell centers and nodes arrays
    X  = spread(xlin, dim=2, ncopies=N)     ! x cell center
    Y  = spread(xlin, dim=1, ncopies=N)     ! y cell center
    
    ! Set initial conditions for desired problem   
    call initialize_problem()

    ! Calculate initial conserved variables
    call get_conserved(rho, vx, vy, P, Bx, By, gamma, vol, N, N, Mass, Momx, Momy, Energy)

    ! Cleanup
    if (allocated(X)) deallocate(X)
    if (allocated(Y)) deallocate(Y)
    if (allocated(xlin)) deallocate(xlin)
    if (allocated(Az)) deallocate(Az)
    
    
    !-----------------
    !--- Main loop ---
    !-----------------
    t = 0.d0           ! Init sim time
    outputCount = 1    ! Init output counter 
    timeStampCount = 1 ! Init timestamp counter
    do while (t < tEnd)

        ! Update primitive vars from conserved vars and face-centered B-fields
        call average_face_to_cell_B(b_x, b_y, N, N, Bx, By)
        call get_primitive(Mass, Momx, Momy, Energy, Bx, By, &
                           gamma, vol, N, N, rho, vx, vy, P)
        inv_rho = 1.0d0 / rho

        ! Calc local wave speeds for CFL-limited timestep:
        c0_sq = gamma * (P - 0.5d0 * (Bx*Bx + By*By)) * inv_rho       ! (sound speed)^2
        ca_sq = (Bx*Bx + By*By) * inv_rho                             ! (Alfven speed)^2
        cf    = sqrt(c0_sq + ca_sq)                                   ! fast mag.sonic speed
        where (rho < 1.0d-3)
            cf = min(cf, 50.0d0)                                      ! Set bound on cf 
        end where

        ! Calculate timestep from wave speeds
        dt    = courant_fac * minval(dx / (cf + sqrt(vx*vx + vy*vy))) ! timestep

        ! Reduce dt if approaching next output time
        if (t + dt > outputCount * tOut) then
            dt = outputCount * tOut - t
        end if
        dt = max(dt, 1e-12)
        
        ! Fill ghost cells for all primitives, then compute gradients.
        ! fill_ghost_cells pads each N x N field to (N+2) x (N+2) using
        ! BC flags from mhd_config. The padded arrays are used by both
        ! compute_gradients and apply_slope_limiter below.
        call fill_ghost_cells(rho, rho_pad)
        call fill_ghost_cells(vx,  vx_pad)
        call fill_ghost_cells(vy,  vy_pad)
        call fill_ghost_cells(P,   P_pad)
        call fill_ghost_cells(Bx,  Bx_pad)
        call fill_ghost_cells(By,  By_pad)

        call compute_gradients(rho_pad, dx, N, N, rho_dx, rho_dy)
        call compute_gradients(vx_pad,  dx, N, N, vx_dx,  vx_dy)
        call compute_gradients(vy_pad,  dx, N, N, vy_dx,  vy_dy)
        call compute_gradients(P_pad,   dx, N, N, P_dx,   P_dy)
        call compute_gradients(Bx_pad,  dx, N, N, Bx_dx,  Bx_dy)
        call compute_gradients(By_pad,  dx, N, N, By_dx,  By_dy)

        ! Apply slope limiter to gradients if enabled. Helps limit
        ! oscillations near discontinuities
        if (useSlopeLimiting) then
            call apply_slope_limiter(rho_pad, dx, N, N, rho_dx, rho_dy)
            call apply_slope_limiter(vx_pad,  dx, N, N, vx_dx,  vx_dy)
            call apply_slope_limiter(vy_pad,  dx, N, N, vy_dx,  vy_dy)
            call apply_slope_limiter(P_pad,   dx, N, N, P_dx,   P_dy)
            call apply_slope_limiter(Bx_pad,  dx, N, N, Bx_dx,  Bx_dy)
            call apply_slope_limiter(By_pad,  dx, N, N, By_dx,  By_dy)
        end if

        ! Extrapolate half-step in time (prediction step)
        !   - Advance primitive vars by dt/2
        !   - Use gradients and velocities to predict vals @ time t + dt/2
        !   - Predicted values (rho_prime, vx_prime, etc.) are used for
        !     boundary reconstruction and flux calculations

        ! Density prediction (w/ advection and compression)
        rho_prime = rho - 0.5d0 * dt * (vx*rho_dx + rho*vx_dx + vy*rho_dy + rho*vy_dy)
        rho_prime = max(rho_prime, rho_floor)

        ! X velocity prediction (w/ advection, pressure grad, and Lorentz force)
        vx_prime = vx - 0.5d0 * dt * ( vx * vx_dx + vy * vx_dy + inv_rho * P_dx &
                      - (2.d0 * Bx * inv_rho) * Bx_dx - (By * inv_rho) * Bx_dy  &
                      - (Bx * inv_rho) * By_dy )

        ! Y velocity prediction (w/ advection, pressure grad, and Lorentz force)
        vy_prime = vy - 0.5d0 * dt * ( vx * vy_dx + vy * vy_dy + inv_rho * P_dy &
                      - (2.d0 * By * inv_rho) * By_dy - (Bx * inv_rho) * By_dx  &
                      - (By * inv_rho) * Bx_dx )

        ! Pressure prediction (w/ advection, compresstion, and magnetic contributions)
        P_prime = P - 0.5d0 * dt * ((gamma * (P - 0.5d0*(Bx*Bx + By*By)) + By*By)*vx_dx &
                    - Bx*By*vy_dx + vx*P_dx + (gamma-2.d0) * (Bx*vx + By*vy) * Bx_dx    &
                    - By*Bx*vx_dy + (gamma * (P - 0.5d0*(Bx*Bx + By*By)) + Bx*Bx)*vy_dy &
                    + vy * P_dy + (gamma - 2.d0) * (Bx * vx + By * vy) * By_dy )
        ! NOTE: P_prime floor is applied after the B predictions below so that
        ! Bx_prime/By_prime can be used. Flooring before would use full-step B,
        ! which is inconsistent with the half-step total pressure formulation.

        ! Bx prediction (induction equation for Bx)
        Bx_prime = Bx - 0.5d0 * dt * (-By * vx_dy + Bx * vy_dy + vy * Bx_dy - vx * By_dy)

        ! By prediction (induction equation for By)
        By_prime = By - 0.5d0 * dt * ( By * vx_dx - Bx * vy_dx - vy * Bx_dx + vx * By_dx)

        ! floor P_prime (total pressure) via the predicted B fields to ensure pos. thermal P
        P_prime = max(P_prime, 0.5d0*(Bx_prime*Bx_prime + By_prime*By_prime) + P_floor)


        ! Extrapolate in space to face centers
        !   - Uses the predicted vars and their grads to reconstruct left/right states at 
        !     each cell face in both x and y directions (MUSCL reconstruction scheme)
        call reconstruction(rho_prime,rho_dx,rho_dy,dx,N,N,rho_XL,rho_XR,rho_YL,rho_YR, .true.)
        call reconstruction(vx_prime, vx_dx, vx_dy, dx,N,N,vx_XL, vx_XR, vx_YL, vx_YR, .false.)
        call reconstruction(vy_prime, vy_dx, vy_dy, dx,N,N,vy_XL, vy_XR, vy_YL, vy_YR, .false.)
        call reconstruction(P_prime,  P_dx,  P_dy,  dx,N,N,P_XL,  P_XR,  P_YL,  P_YR,   .true.)
        call reconstruction(Bx_prime, Bx_dx, Bx_dy, dx,N,N,Bx_XL, Bx_XR, Bx_YL, Bx_YR, .false.)
        call reconstruction(By_prime, By_dx, By_dy, dx,N,N,By_XL, By_XR, By_YL, By_YR, .false.)

        ! Thermal pressure positivity check (MOOD extension):
        ! Fall back to cell-centre on any face where P - 0.5*B^2 < P_floor
        call thermal_pressure_check(N, N,                              &
            P_XL,  P_XR,  P_YL,  P_YR,                               &
            Bx_XL, Bx_XR, Bx_YL, Bx_YR,                             &
            By_XL, By_XR, By_YL, By_YR,                              &
            rho_XL, rho_XR, rho_YL, rho_YR,                          &
            vx_XL,  vx_XR,  vx_YL,  vx_YR,                          &
            vy_XL,  vy_XR,  vy_YL,  vy_YR,                          &
            rho_prime, vx_prime, vy_prime, P_prime, Bx_prime, By_prime)


        ! Compute fluxes of conserved vars (local Lax-Friedrichs/Rusanov)
        !   - For x-faces: use left/right states in x-direction
        !   - For y-faces: use left/right states in y-direction
        call compute_fluxes(rho_XL, rho_XR, vx_XL, vx_XR, vy_XL, vy_XR, P_XL, P_XR, &
                            Bx_XL, Bx_XR, By_XL, By_XR, gamma, N, N, flux_Mass_X,   &
                            flux_Momx_X, flux_Momy_X, flux_Energy_X, flux_By_X)
        call compute_fluxes(rho_YL, rho_YR, vy_YL, vy_YR, vx_YL, vx_YR, P_YL, P_YR, &
                            By_YL, By_YR, Bx_YL, Bx_YR, gamma, N, N, flux_Mass_Y,   &
                            flux_Momy_Y, flux_Momx_Y, flux_Energy_Y, flux_Bx_Y)


        ! Update conserved variables using computed fluxes and apply 
        ! constrained transport update for face-centered B-fields
        call update_conserved(Mass, flux_Mass_X, flux_Mass_Y, dx, dt, N, N)
        call update_conserved(Momx, flux_Momx_X, flux_Momx_Y, dx, dt, N, N)
        call update_conserved(Momy, flux_Momy_X, flux_Momy_Y, dx, dt, N, N)
        call update_conserved(Energy, flux_Energy_X, flux_Energy_Y, dx, dt, N, N)
        call apply_conserved_floors(Mass, Momx, Momy, Energy, Bx, By, vol, N, N)
        call constrained_transport(b_x, b_y, flux_By_X, flux_Bx_Y, dx, dt, N, N)

        ! Outflow boundary flux correction for xlo and ylo.
        ! update_conserved applies interior face fluxes but adds no incoming flux
        ! from outside at non-periodic low-side boundaries. This block adds the
        ! one-sided MHD flux at those faces
        if (bc_xlo /= BC_PERIODIC) then
            block
                real(8) :: halfB2_1(N), en_1(N)
                halfB2_1 = 0.5d0*(Bx_prime(1,:)**2 + By_prime(1,:)**2)
                en_1 = (P_prime(1,:) - halfB2_1)/(gamma - 1.d0) &
                     + 0.5d0*rho_prime(1,:)*(vx_prime(1,:)**2 + vy_prime(1,:)**2) + halfB2_1
                Mass(1,:)   = Mass(1,:)   + dt*dx * rho_prime(1,:)*vx_prime(1,:)
                Momx(1,:)   = Momx(1,:)   + dt*dx * (rho_prime(1,:)*vx_prime(1,:)**2 &
                              + P_prime(1,:) - Bx_prime(1,:)**2)
                Momy(1,:)   = Momy(1,:)   + dt*dx * (rho_prime(1,:)*vx_prime(1,:)*vy_prime(1,:) &
                              - Bx_prime(1,:)*By_prime(1,:))
                Energy(1,:) = Energy(1,:) + dt*dx * ((en_1 + P_prime(1,:))*vx_prime(1,:) &
                              - Bx_prime(1,:)*(Bx_prime(1,:)*vx_prime(1,:) + By_prime(1,:)*vy_prime(1,:)))
            end block
        end if

        if (bc_ylo /= BC_PERIODIC) then
            block
                real(8) :: halfB2_1(N), en_1(N)
                halfB2_1 = 0.5d0*(Bx_prime(:,1)**2 + By_prime(:,1)**2)
                en_1 = (P_prime(:,1) - halfB2_1)/(gamma - 1.d0) &
                     + 0.5d0*rho_prime(:,1)*(vx_prime(:,1)**2 + vy_prime(:,1)**2) + halfB2_1
                Mass(:,1)   = Mass(:,1)   + dt*dx * rho_prime(:,1)*vy_prime(:,1)
                Momy(:,1)   = Momy(:,1)   + dt*dx * (rho_prime(:,1)*vy_prime(:,1)**2 &
                              + P_prime(:,1) - By_prime(:,1)**2)
                Momx(:,1)   = Momx(:,1)   + dt*dx * (rho_prime(:,1)*vy_prime(:,1)*vx_prime(:,1) &
                              - By_prime(:,1)*Bx_prime(:,1))
                Energy(:,1) = Energy(:,1) + dt*dx * ((en_1 + P_prime(:,1))*vy_prime(:,1) &
                              - By_prime(:,1)*(By_prime(:,1)*vy_prime(:,1) + Bx_prime(:,1)*vx_prime(:,1)))
            end block
        end if

        ! Update time
        t = t + dt

        ! Calc div B to see if it's near 0 (should be)
        call compute_divB(b_x, b_y, dx, N, N, divB)

        ! Output data every interval specified in 
        ! config file, and again at the end of the sim
        if (t >= outputCount * tOut .or. t >= tEnd) then
            
            ! Update primitives 
            call average_face_to_cell_B(b_x, b_y, N, N, Bx, By)
            call get_primitive(Mass, Momx, Momy, Energy, Bx, By, &
                               gamma, vol, N, N, rho, vx, vy, P)

            ! Append the current field data
            rho_all(outputCount, :, :) = rho
            P_all(outputCount, :, :) = P - 0.5d0*(Bx*Bx + By*By)
            Bx_all(outputCount, :, :) = Bx
            By_all(outputCount, :, :) = By
            Vx_all(outputCount, :, :) = vx
            Vy_all(outputCount, :, :) = vy
            time_all(outputCount) = t
            outputCount = outputCount + 1
        end if

        ! Print sim progress and average B-field divergence to terminal.
        ! Skip this if doing Monte Carlo
        if (t >= timeStampCount * outInterval .and. t < tEnd + 1.0d-8 .and. problem_type .ne. 5) then
            print '(I3, A, 2X, A, F7.3, 2X, A, ES12.4)', &
                NINT(100.0d0 * (timeStampCount * outInterval) / tEnd), '%:', &
                't =', t, ', avg(|divB|) =', sum(abs(divB)) / dble(N*N)
            timeStampCount = timeStampCount + 1
        end if   

    end do

    ! Trim arrays
    actual_nsnap = outputCount - 1 
    allocate(rho_trim(actual_nsnap, N, N))
    allocate(P_trim(actual_nsnap, N, N))
    allocate(Bx_trim(actual_nsnap, N, N))
    allocate(By_trim(actual_nsnap, N, N))
    allocate(Vx_trim(actual_nsnap, N, N))
    allocate(Vy_trim(actual_nsnap, N, N))
    allocate(time_trim(actual_nsnap))
    rho_trim = rho_all(1:actual_nsnap, :, :)
    P_trim   = P_all(1:actual_nsnap, :, :)
    Bx_trim  = Bx_all(1:actual_nsnap, :, :)
    By_trim  = By_all(1:actual_nsnap, :, :)
    Vx_trim  = Vx_all(1:actual_nsnap, :, :)
    Vy_trim  = Vy_all(1:actual_nsnap, :, :)
    time_trim = time_all(1:actual_nsnap)

    ! Write trimmed arrays
    full_h5_path = trim(out_path) // trim(h5_filename)
    call write_prims_to_hdf5(full_h5_path, rho_trim, P_trim, Bx_trim, By_trim, & 
                            Vx_trim, Vy_trim, time_trim, gamma_actual,         &
                            M_s_actual, beta_actual)

    ! Clean up
    deallocate(rho_trim, P_trim, Bx_trim, By_trim, Vx_trim, Vy_trim, time_trim)

end program main
