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

    ! For trimmed output tensors
    integer :: actual_nsnap
    real(kind=8), dimension(:,:,:), allocatable :: rho_trim, P_trim, Bx_trim, By_trim, Vx_trim, Vy_trim
    real(kind=8), dimension(:), allocatable :: time_trim
    character(len=512) :: full_h5_path

    ! SSP-RK3 stage storage: conserved state at the start of each timestep (u^n).
    real(kind=8), allocatable :: Mass_0(:,:), Momx_0(:,:), Momy_0(:,:), Energy_0(:,:)
    real(kind=8), allocatable :: bx_0(:,:), by_0(:,:)

    ! Get desired problem, output path, and filename
    call get_user_options()

    ! Setup simulation
    call initialize_variables()

    ! Allocate RK3 stage storage
    allocate(Mass_0(N,N), Momx_0(N,N), Momy_0(N,N), Energy_0(N,N))
    allocate(bx_0(N,N), by_0(N,N))

    ! Generate uniform grid
    dx   = boxsize / dble(N)                ! cell width
    vol  = dx*dx                            ! cell area
    xlin = [( (i - 0.5d0), i = 1, N)] * dx  ! cell centers

    ! Fill 2D grid cell centers and nodes arrays
    X  = spread(xlin, dim=2, ncopies=N)     ! x cell center
    Y  = spread(xlin, dim=1, ncopies=N)     ! y cell center

    ! Set initial conditions for desired problem
    call initialize_problem()

    ! Ensure gamma is written to HDF5
    gamma_actual = gamma

    ! Calculate initial conserved variables
    call get_conserved(rho, vx, vy, P, Bx, By, gamma, vol, N, N, Mass, Momx, Momy, Energy)

    ! Cleanup
    if (allocated(X))    deallocate(X)
    if (allocated(Y))    deallocate(Y)
    if (allocated(xlin)) deallocate(xlin)
    if (allocated(Az))   deallocate(Az)


    !-----------------
    !--- Main loop ---
    !-----------------
    t              = 0.d0
    outputCount    = 1
    timeStampCount = 1

    do while (t < tEnd)

        ! Compute primitive variables from current conserved state.
        ! Bx, By are updated here from staggered face fields b_x, b_y.
        call average_face_to_cell_B(b_x, b_y, N, N, Bx, By)
        call get_primitive(Mass, Momx, Momy, Energy, Bx, By, &
                           gamma, vol, N, N, rho, vx, vy, P)
        inv_rho = 1.0d0 / rho

        ! CFL-limited timestep (calc once per timestep from u^n)
        c0_sq = gamma * (P - 0.5d0 * (Bx*Bx + By*By)) * inv_rho       ! sound speed^2
        ca_sq = (Bx*Bx + By*By) * inv_rho                             ! Alfven speed^2
        cf    = sqrt(c0_sq + ca_sq)                                   ! fast magnetosonic speed
        where (rho < 1.0d-3) cf = min(cf, 50.0d0)                     ! CFL guard for near-vacuum
        dt    = courant_fac * minval(dx / (cf + sqrt(vx*vx + vy*vy)))

        ! Reduce dt if approaching next output time
        if (t + dt > outputCount * tOut) dt = outputCount * tOut - t
        dt = max(dt, 1e-12)

        ! Save u^n for SSP-RK3 linear combinations.
        Mass_0   = Mass;   Momx_0 = Momx; Momy_0 = Momy
        Energy_0 = Energy; bx_0   = b_x;  by_0   = b_y

        !================================================================
        ! SSP-RK3 time integration  (Shu & Osher 1988):
        !
        !   Stage 1:  u^(1)   = u^n + dt * L(u^n)
        !   Stage 2:  u^(2)   = 3/4 u^n + 1/4 [ u^(1) + dt * L(u^(1)) ]
        !   Stage 3:  u^(n+1) = 1/3 u^n + 2/3 [ u^(2) + dt * L(u^(2)) ]
        !
        ! L(u) is spatial operator: ghost cells -> gradients -> slope
        ! limiter -> reconstruction -> thermal pressure check -> fluxes.
        ! Implemented in subroutine compute_stage_fluxes
        !================================================================


        ! ============================================================
        ! Stage 1:  u^(1) = u^n + dt * L(u^n)
        ! Primitives rho/vx/vy/P/Bx/By already reflect u^n from the
        ! get_primitive call above - no extra work needed here.
        ! ============================================================
        call compute_stage_fluxes()

        call update_conserved(Mass,   flux_Mass_X,   flux_Mass_Y,   dx, dt, N, N)
        call update_conserved(Momx,   flux_Momx_X,   flux_Momx_Y,   dx, dt, N, N)
        call update_conserved(Momy,   flux_Momy_X,   flux_Momy_Y,   dx, dt, N, N)
        call update_conserved(Energy, flux_Energy_X, flux_Energy_Y, dx, dt, N, N)
        call apply_conserved_floors(Mass, Momx, Momy, Energy, Bx, By, vol, gamma, N, N)
        call constrained_transport(b_x, b_y, flux_By_X, flux_Bx_Y, dx, dt, N, N)

        ! Outflow boundary flux correction (xlo and ylo non-periodic faces).
        ! Uses cell-center primitives at the current stage (rho, vx, vy, P, Bx, By).
        if (bc_xlo /= BC_PERIODIC) then
            block
                real(8) :: halfB2_1(N), en_1(N)
                halfB2_1 = 0.5d0*(Bx(1,:)**2 + By(1,:)**2)
                en_1 = (P(1,:) - halfB2_1)/(gamma - 1.d0) &
                     + 0.5d0*rho(1,:)*(vx(1,:)**2 + vy(1,:)**2) + halfB2_1
                Mass(1,:)   = Mass(1,:)   + dt*dx * rho(1,:)*vx(1,:)
                Momx(1,:)   = Momx(1,:)   + dt*dx * (rho(1,:)*vx(1,:)**2 &
                              + P(1,:) - Bx(1,:)**2)
                Momy(1,:)   = Momy(1,:)   + dt*dx * (rho(1,:)*vx(1,:)*vy(1,:) &
                              - Bx(1,:)*By(1,:))
                Energy(1,:) = Energy(1,:) + dt*dx * ((en_1 + P(1,:))*vx(1,:) &
                              - Bx(1,:)*(Bx(1,:)*vx(1,:) + By(1,:)*vy(1,:)))
            end block
        end if

        if (bc_ylo /= BC_PERIODIC) then
            block
                real(8) :: halfB2_1(N), en_1(N)
                halfB2_1 = 0.5d0*(Bx(:,1)**2 + By(:,1)**2)
                en_1 = (P(:,1) - halfB2_1)/(gamma - 1.d0) &
                     + 0.5d0*rho(:,1)*(vx(:,1)**2 + vy(:,1)**2) + halfB2_1
                Mass(:,1)   = Mass(:,1)   + dt*dx * rho(:,1)*vy(:,1)
                Momy(:,1)   = Momy(:,1)   + dt*dx * (rho(:,1)*vy(:,1)**2 &
                              + P(:,1) - By(:,1)**2)
                Momx(:,1)   = Momx(:,1)   + dt*dx * (rho(:,1)*vy(:,1)*vx(:,1) &
                              - By(:,1)*Bx(:,1))
                Energy(:,1) = Energy(:,1) + dt*dx * ((en_1 + P(:,1))*vy(:,1) &
                              - By(:,1)*(By(:,1)*vy(:,1) + Bx(:,1)*vx(:,1)))
            end block
        end if
        ! Conserved vars now hold u^(1)


        ! ============================================================
        ! Stage 2:  u^(2) = 3/4 u^n  +  1/4 [ u^(1) + dt * L(u^(1)) ]
        ! ============================================================
        call average_face_to_cell_B(b_x, b_y, N, N, Bx, By)
        call get_primitive(Mass, Momx, Momy, Energy, Bx, By, &
                           gamma, vol, N, N, rho, vx, vy, P)
        inv_rho = 1.0d0 / rho

        call compute_stage_fluxes()

        call update_conserved(Mass,   flux_Mass_X,   flux_Mass_Y,   dx, dt, N, N)
        call update_conserved(Momx,   flux_Momx_X,   flux_Momx_Y,   dx, dt, N, N)
        call update_conserved(Momy,   flux_Momy_X,   flux_Momy_Y,   dx, dt, N, N)
        call update_conserved(Energy, flux_Energy_X, flux_Energy_Y, dx, dt, N, N)
        call apply_conserved_floors(Mass, Momx, Momy, Energy, Bx, By, vol, gamma, N, N)
        call constrained_transport(b_x, b_y, flux_By_X, flux_Bx_Y, dx, dt, N, N)

        if (bc_xlo /= BC_PERIODIC) then
            block
                real(8) :: halfB2_1(N), en_1(N)
                halfB2_1 = 0.5d0*(Bx(1,:)**2 + By(1,:)**2)
                en_1 = (P(1,:) - halfB2_1)/(gamma - 1.d0) &
                     + 0.5d0*rho(1,:)*(vx(1,:)**2 + vy(1,:)**2) + halfB2_1
                Mass(1,:)   = Mass(1,:)   + dt*dx * rho(1,:)*vx(1,:)
                Momx(1,:)   = Momx(1,:)   + dt*dx * (rho(1,:)*vx(1,:)**2 &
                              + P(1,:) - Bx(1,:)**2)
                Momy(1,:)   = Momy(1,:)   + dt*dx * (rho(1,:)*vx(1,:)*vy(1,:) &
                              - Bx(1,:)*By(1,:))
                Energy(1,:) = Energy(1,:) + dt*dx * ((en_1 + P(1,:))*vx(1,:) &
                              - Bx(1,:)*(Bx(1,:)*vx(1,:) + By(1,:)*vy(1,:)))
            end block
        end if

        if (bc_ylo /= BC_PERIODIC) then
            block
                real(8) :: halfB2_1(N), en_1(N)
                halfB2_1 = 0.5d0*(Bx(:,1)**2 + By(:,1)**2)
                en_1 = (P(:,1) - halfB2_1)/(gamma - 1.d0) &
                     + 0.5d0*rho(:,1)*(vx(:,1)**2 + vy(:,1)**2) + halfB2_1
                Mass(:,1)   = Mass(:,1)   + dt*dx * rho(:,1)*vy(:,1)
                Momy(:,1)   = Momy(:,1)   + dt*dx * (rho(:,1)*vy(:,1)**2 &
                              + P(:,1) - By(:,1)**2)
                Momx(:,1)   = Momx(:,1)   + dt*dx * (rho(:,1)*vy(:,1)*vx(:,1) &
                              - By(:,1)*Bx(:,1))
                Energy(:,1) = Energy(:,1) + dt*dx * ((en_1 + P(:,1))*vy(:,1) &
                              - By(:,1)*(By(:,1)*vy(:,1) + Bx(:,1)*vx(:,1)))
            end block
        end if

        ! Blend:  u^(2) = 3/4 u^n  +  1/4 (u^(1) + dt*L(u^(1)))
        Mass   = 0.75d0*Mass_0   + 0.25d0*Mass
        Momx   = 0.75d0*Momx_0   + 0.25d0*Momx
        Momy   = 0.75d0*Momy_0   + 0.25d0*Momy
        Energy = 0.75d0*Energy_0 + 0.25d0*Energy
        b_x    = 0.75d0*bx_0     + 0.25d0*b_x
        b_y    = 0.75d0*by_0     + 0.25d0*b_y
        call apply_conserved_floors(Mass, Momx, Momy, Energy, Bx, By, vol, gamma, N, N)


        ! ============================================================
        ! Stage 3:  u^(n+1) = 1/3 u^n  +  2/3 [ u^(2) + dt * L(u^(2)) ]
        ! ============================================================
        call average_face_to_cell_B(b_x, b_y, N, N, Bx, By)
        call get_primitive(Mass, Momx, Momy, Energy, Bx, By, &
                           gamma, vol, N, N, rho, vx, vy, P)
        inv_rho = 1.0d0 / rho

        call compute_stage_fluxes()

        call update_conserved(Mass,   flux_Mass_X,   flux_Mass_Y,   dx, dt, N, N)
        call update_conserved(Momx,   flux_Momx_X,   flux_Momx_Y,   dx, dt, N, N)
        call update_conserved(Momy,   flux_Momy_X,   flux_Momy_Y,   dx, dt, N, N)
        call update_conserved(Energy, flux_Energy_X, flux_Energy_Y, dx, dt, N, N)
        call apply_conserved_floors(Mass, Momx, Momy, Energy, Bx, By, vol, gamma, N, N)
        call constrained_transport(b_x, b_y, flux_By_X, flux_Bx_Y, dx, dt, N, N)

        if (bc_xlo /= BC_PERIODIC) then
            block
                real(8) :: halfB2_1(N), en_1(N)
                halfB2_1 = 0.5d0*(Bx(1,:)**2 + By(1,:)**2)
                en_1 = (P(1,:) - halfB2_1)/(gamma - 1.d0) &
                     + 0.5d0*rho(1,:)*(vx(1,:)**2 + vy(1,:)**2) + halfB2_1
                Mass(1,:)   = Mass(1,:)   + dt*dx * rho(1,:)*vx(1,:)
                Momx(1,:)   = Momx(1,:)   + dt*dx * (rho(1,:)*vx(1,:)**2 &
                              + P(1,:) - Bx(1,:)**2)
                Momy(1,:)   = Momy(1,:)   + dt*dx * (rho(1,:)*vx(1,:)*vy(1,:) &
                              - Bx(1,:)*By(1,:))
                Energy(1,:) = Energy(1,:) + dt*dx * ((en_1 + P(1,:))*vx(1,:) &
                              - Bx(1,:)*(Bx(1,:)*vx(1,:) + By(1,:)*vy(1,:)))
            end block
        end if

        if (bc_ylo /= BC_PERIODIC) then
            block
                real(8) :: halfB2_1(N), en_1(N)
                halfB2_1 = 0.5d0*(Bx(:,1)**2 + By(:,1)**2)
                en_1 = (P(:,1) - halfB2_1)/(gamma - 1.d0) &
                     + 0.5d0*rho(:,1)*(vx(:,1)**2 + vy(:,1)**2) + halfB2_1
                Mass(:,1)   = Mass(:,1)   + dt*dx * rho(:,1)*vy(:,1)
                Momy(:,1)   = Momy(:,1)   + dt*dx * (rho(:,1)*vy(:,1)**2 &
                              + P(:,1) - By(:,1)**2)
                Momx(:,1)   = Momx(:,1)   + dt*dx * (rho(:,1)*vy(:,1)*vx(:,1) &
                              - By(:,1)*Bx(:,1))
                Energy(:,1) = Energy(:,1) + dt*dx * ((en_1 + P(:,1))*vy(:,1) &
                              - By(:,1)*(By(:,1)*vy(:,1) + Bx(:,1)*vx(:,1)))
            end block
        end if

        ! Blend:  u^(n+1) = 1/3 u^n  +  2/3 (u^(2) + dt*L(u^(2)))
        Mass   = (1.d0/3.d0)*Mass_0   + (2.d0/3.d0)*Mass
        Momx   = (1.d0/3.d0)*Momx_0   + (2.d0/3.d0)*Momx
        Momy   = (1.d0/3.d0)*Momy_0   + (2.d0/3.d0)*Momy
        Energy = (1.d0/3.d0)*Energy_0 + (2.d0/3.d0)*Energy
        b_x    = (1.d0/3.d0)*bx_0     + (2.d0/3.d0)*b_x
        b_y    = (1.d0/3.d0)*by_0     + (2.d0/3.d0)*b_y
        call apply_conserved_floors(Mass, Momx, Momy, Energy, Bx, By, vol, gamma, N, N)


        ! ----------------------------------------------------------------
        ! Advance time, diagnostics, and output
        ! ----------------------------------------------------------------
        t = t + dt

        ! Calc div B to see if it's near 0 (should be)
        call compute_divB(b_x, b_y, dx, N, N, divB)

        ! Output data every interval specified in config file, and again at the end
        if (t >= outputCount * tOut .or. t >= tEnd) then
            call average_face_to_cell_B(b_x, b_y, N, N, Bx, By)
            call get_primitive(Mass, Momx, Momy, Energy, Bx, By, &
                               gamma, vol, N, N, rho, vx, vy, P)

            rho_all(outputCount, :, :) = rho
            P_all(outputCount, :, :)   = P - 0.5d0*(Bx*Bx + By*By)
            Bx_all(outputCount, :, :)  = Bx
            By_all(outputCount, :, :)  = By
            Vx_all(outputCount, :, :)  = vx
            Vy_all(outputCount, :, :)  = vy
            time_all(outputCount)      = t
            outputCount = outputCount + 1
        end if

        ! Print sim progress and average B-field divergence.
        ! Skip this if doing Monte Carlo.
        if (t >= timeStampCount * outInterval .and. t < tEnd + 1.0d-8 .and. problem_type .ne. 5) then
            print '(I3, A, 2X, A, F7.3, 2X, A, ES12.4)', &
                NINT(100.0d0 * (timeStampCount * outInterval) / tEnd), '%:', &
                't =', t, ', avg(|divB|) =', sum(abs(divB)) / dble(N*N)
            timeStampCount = timeStampCount + 1
        end if

    end do

    ! ----------------------------------------------------------------
    ! Trim arrays to the actual number of snapshots written and output.
    ! ----------------------------------------------------------------
    actual_nsnap = outputCount - 1
    allocate(rho_trim(actual_nsnap, N, N))
    allocate(P_trim(actual_nsnap, N, N))
    allocate(Bx_trim(actual_nsnap, N, N))
    allocate(By_trim(actual_nsnap, N, N))
    allocate(Vx_trim(actual_nsnap, N, N))
    allocate(Vy_trim(actual_nsnap, N, N))
    allocate(time_trim(actual_nsnap))
    rho_trim  = rho_all(1:actual_nsnap, :, :)
    P_trim    = P_all(1:actual_nsnap, :, :)
    Bx_trim   = Bx_all(1:actual_nsnap, :, :)
    By_trim   = By_all(1:actual_nsnap, :, :)
    Vx_trim   = Vx_all(1:actual_nsnap, :, :)
    Vy_trim   = Vy_all(1:actual_nsnap, :, :)
    time_trim = time_all(1:actual_nsnap)

    full_h5_path = trim(out_path) // trim(h5_filename)
    call write_prims_to_hdf5(full_h5_path, rho_trim, P_trim, Bx_trim, By_trim, &
                             Vx_trim, Vy_trim, time_trim, gamma_actual,         &
                             M_s_actual, beta_actual)

    deallocate(rho_trim, P_trim, Bx_trim, By_trim, Vx_trim, Vy_trim, time_trim)
    deallocate(Mass_0, Momx_0, Momy_0, Energy_0, bx_0, by_0)

contains


    ! ====================================================================
    ! compute_stage_fluxes
    !
    ! Evaluates spatial operator L(u) for one SSP-RK3 stage. Takes in 
    ! primitive arrays (rho, vx, vy, P, Bx, By).
    !
    ! Fills:  flux_Mass_X/Y, flux_Momx_X/Y, flux_Momy_X/Y,
    !         flux_Energy_X/Y, flux_By_X, flux_Bx_Y
    !
    ! Reconstruction uses the cell-center primitives directly
    ! ====================================================================
    subroutine compute_stage_fluxes()

        ! Ghost cells for all six primitives
        call fill_ghost_cells(rho, rho_pad)
        call fill_ghost_cells(vx,  vx_pad)
        call fill_ghost_cells(vy,  vy_pad)
        call fill_ghost_cells(P,   P_pad)
        call fill_ghost_cells(Bx,  Bx_pad)
        call fill_ghost_cells(By,  By_pad)

        ! Second-order centered gradients
        call compute_gradients(rho_pad, dx, N, N, rho_dx, rho_dy)
        call compute_gradients(vx_pad,  dx, N, N, vx_dx,  vx_dy)
        call compute_gradients(vy_pad,  dx, N, N, vy_dx,  vy_dy)
        call compute_gradients(P_pad,   dx, N, N, P_dx,   P_dy)
        call compute_gradients(Bx_pad,  dx, N, N, Bx_dx,  Bx_dy)
        call compute_gradients(By_pad,  dx, N, N, By_dx,  By_dy)

        ! Slope limiter
        if (useSlopeLimiting) then
            call apply_slope_limiter(rho_pad, dx, N, N, rho_dx, rho_dy)
            call apply_slope_limiter(vx_pad,  dx, N, N, vx_dx,  vx_dy)
            call apply_slope_limiter(vy_pad,  dx, N, N, vy_dx,  vy_dy)
            call apply_slope_limiter(P_pad,   dx, N, N, P_dx,   P_dy)
            call apply_slope_limiter(Bx_pad,  dx, N, N, Bx_dx,  Bx_dy)
            call apply_slope_limiter(By_pad,  dx, N, N, By_dx,  By_dy)
        end if

        ! MUSCL reconstruction to cell faces using cell-centre values.
        call reconstruction(rho, rho_dx, rho_dy, dx, N, N, rho_XL, rho_XR, rho_YL, rho_YR,  .true.)
        call reconstruction(vx,  vx_dx,  vx_dy,  dx, N, N, vx_XL,  vx_XR,  vx_YL,  vx_YR,  .false.)
        call reconstruction(vy,  vy_dx,  vy_dy,  dx, N, N, vy_XL,  vy_XR,  vy_YL,  vy_YR,  .false.)
        call reconstruction(P,   P_dx,   P_dy,   dx, N, N, P_XL,   P_XR,   P_YL,   P_YR,    .true.)
        call reconstruction(Bx,  Bx_dx,  Bx_dy,  dx, N, N, Bx_XL,  Bx_XR,  Bx_YL,  Bx_YR,  .false.)
        call reconstruction(By,  By_dx,  By_dy,  dx, N, N, By_XL,  By_XR,  By_YL,  By_YR,  .false.)

        ! Thermal pressure positivity check (MOOD extension).
        ! Cell-centre values are the fallback (was _prime in Hancock scheme).
        call thermal_pressure_check(N, N,   &
            P_XL,  P_XR,  P_YL,  P_YR,      &
            Bx_XL, Bx_XR, Bx_YL, Bx_YR,     &
            By_XL, By_XR, By_YL, By_YR,     &
            rho_XL, rho_XR, rho_YL, rho_YR, &
            vx_XL,  vx_XR,  vx_YL,  vx_YR,  &
            vy_XL,  vy_XR,  vy_YL,  vy_YR,  &
            rho, vx, vy, P, Bx, By)

        ! Riemann fluxes in x-direction (normal = vx, tangential = vy)
        call compute_fluxes(rho_XL, rho_XR, vx_XL, vx_XR, vy_XL, vy_XR, P_XL, P_XR, &
                            Bx_XL, Bx_XR, By_XL, By_XR, gamma, N, N, flux_Mass_X,   &
                            flux_Momx_X, flux_Momy_X, flux_Energy_X, flux_By_X)

        ! Riemann fluxes in y-direction (rotated: vy=normal, vx=tangential, By=normal B)
        call compute_fluxes(rho_YL, rho_YR, vy_YL, vy_YR, vx_YL, vx_YR, P_YL, P_YR, &
                            By_YL, By_YR, Bx_YL, Bx_YR, gamma, N, N, flux_Mass_Y,   &
                            flux_Momy_Y, flux_Momx_Y, flux_Energy_Y, flux_Bx_Y)

    end subroutine compute_stage_fluxes

end program main
