module mhd_flux
    use mhd_field_ops
    use mhd_config
    !---------------------------------------------------------------------------
    ! Purpose: Subroutines for flux computations in MHD solver: apply numerical
    !          fluxes to conserved fields, constrained transport for divergence
    !          control, and flux function evaluation. Assumes periodic BCs
    !---------------------------------------------------------------------------


    implicit none
    private
    public update_conserved, constrained_transport, compute_fluxes
    
contains


    subroutine update_conserved(F, flux_F_X, flux_F_Y, dx, dt, nx, ny)
    ! 
    !   Apply fluxes to conserved field using periodic BCs. Updates conserved 
    !   variable F by applying the net flux diffs. in both x and y directions
    !
    !   Inputs:
    !       - flux_F_X : x-direction flux (nx, ny)
    !       - flux_F_Y : y-direction flux (nx, ny)
    !       - dx       : cell size
    !       - dt       : time step
    !       - nx       : x grid size
    !       - ny       : y grid size
    !
    !   InOuts:
    !       - F : conserved variable (nx, ny)
    ! 
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: dx, dt
        real(8), intent(in) :: flux_F_X(nx, ny), flux_F_Y(nx, ny)
        real(8), intent(inout) :: F(nx, ny)
        real(8) :: dtdx

        ! Update F in each cell by adding net in/out fluxes
        ! For each cell (i, j):
        !   - subtract outgoing fluxes
        !   - add incoming fluxes from left (i-1) and bottom (j-1)
        dtdx = dt*dx
        F = F - dtdx*flux_F_X + dtdx*cshift(flux_F_X, -1, 1) &
              - dtdx*flux_F_Y + dtdx*cshift(flux_F_Y, -1, 2)

    end subroutine update_conserved


    subroutine constrained_transport(bx, by, flux_By_X, flux_Bx_Y, dx, dt, nx, ny)
    ! 
    !   Apply constrained transport to face-centered magnetic fields, i.e., 
    !   update magnetic field components (bx, by) using discrete curl of the the 
    !   electric field (Ez). Ez is built from B-field fluxes make sure div B ~ 0
    !
    !   Inputs:
    !       - flux_By_X : x-dir flux of By (nx, ny)
    !       - flux_Bx_Y : y-dir flux of Bx (nx, ny)
    !       - dx        : cell size
    !       - dt        : timestep
    !       - nx        : x grid size
    !       - ny        : y grid size
    !
    !   InOuts:
    !       - bx : magnetic field face x component 
    !       - by : magnetic field face y component
    !
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: dx, dt
        real(8), intent(in) :: flux_By_X(nx, ny), flux_Bx_Y(nx, ny)
        real(8), intent(inout) :: bx(nx, ny), by(nx, ny)

        real(8) :: Ez(nx, ny)
        real(8) :: dbx(nx, ny), dby(nx, ny)

        ! Calc z-component of electric field at cell corners
        Ez = 0.25d0 * ( &
            - flux_By_X &
            - cshift(flux_By_X, 1, 2) &   ! (i, j+1) periodic in y
            + flux_Bx_Y &
            + cshift(flux_Bx_Y, 1, 1) )   ! (i+1, j) periodic in x

        ! Calc discrete curl of -Ez and update bx, by
        call compute_curl_2d(-Ez, dx, nx, ny, dbx, dby)
        bx = bx + dt * dbx  ! Update x-face B field
        by = by + dt * dby  ! Update y-face B field

    end subroutine constrained_transport


    subroutine compute_fluxes(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R,   &
                              Bx_L, Bx_R, By_L, By_R, gamma, nx, ny,            &
                              flux_Mass, flux_Momx, flux_Momy, flux_Energy, flux_By)
    !   Compute numerical fluxes for mass, momentum, energy, and transverse magnetic 
    !   field (By) across cell faces. Uses the local Lax-Friedrichs Rusanov flux. 
    !   Basically this averages left and right states and adds a diffusive term 
    !   that's proportional to the maximum local wavespeed to make it stable
    !
    !   Inputs:
    !       - rho_L, rho_R : Left and right densities       [nx, ny]
    !       - vx_L, vx_R   : Left and right x-velocities    [nx, ny]
    !       - vy_L, vy_R   : Left and right y-velocities    [nx, ny]
    !       - P_L, P_R     : Left and right total pressures [nx, ny]
    !       - Bx_L, Bx_R   : Left and right Bx fields       [nx, ny]
    !       - By_L, By_R   : Left and right By fields       [nx, ny]
    !       - gamma        : Adiabatic index (ideal gas)
    !       - nx, ny       : Grid dimensions in x and y
    !
    !   Outputs:
    !       - flux_Mass   : Flux of mass across faces             [nx, ny]
    !       - flux_Momx   : Flux of x-momentum                    [nx, ny]
    !       - flux_Momy   : Flux of y-momentum                    [nx, ny]
    !       - flux_Energy : Flux of energy                        [nx, ny]
    !       - flux_By     : Flux of transverse magnetic field, By [nx, ny]
    !
    !  Method:
    !       1. Compute energies of left and right states.
    !       2. Compute averaged states of primitive and conserved quantities.
    !       3. Evaluate fluxes from the star states.
    !       4. Estimate local maximum wavespeed using fast magnetosonic speed.
    !       5. Apply a stabilizing diffusive term proportional to the wavespeed.
    !
    !  Notes:
    !       - Assumes ideal MHD equations in 2D.
    !       - Uses local Rusanov (LF) flux.
    !       - Magnetic field is split into Bx and By components.
    !
        integer, intent(in) :: nx, ny
        real(8), intent(inout) :: rho_L(nx, ny), rho_R(nx, ny)
        real(8), intent(in) :: vx_L(nx, ny), vx_R(nx, ny)
        real(8), intent(in) :: vy_L(nx, ny), vy_R(nx, ny)
        real(8), intent(inout) :: P_L(nx, ny), P_R(nx, ny)
        real(8), intent(in) :: Bx_L(nx, ny), Bx_R(nx, ny)
        real(8), intent(in) :: By_L(nx, ny), By_R(nx, ny)
        real(8), intent(in) :: gamma

        real(8), intent(out) :: flux_Mass(nx, ny), flux_Momx(nx, ny), flux_Momy(nx, ny)
        real(8), intent(out) :: flux_Energy(nx, ny), flux_By(nx, ny)

        ! Locals arrays
        real(8) :: en_L(nx, ny), en_R(nx, ny), halfBL2(nx, ny), halfBR2(nx, ny)
        real(8) :: rho_avg(nx, ny), inv_rho_avg(nx, ny)
        real(8) :: momx_avg(nx, ny), momy_avg(nx, ny), en_avg(nx, ny)
        real(8) :: Bx_avg(nx, ny), By_avg(nx, ny), P_avg(nx, ny)
        real(8) :: c_L2(nx, ny), c_R2(nx, ny), C_L(nx, ny), C_R(nx, ny), C(nx, ny)

        ! Compute total energy for left and right states
        !   - en = internal + kinetic + magnetic energy
        !       - Internal: (P - 0.5 * B^2) / (gamma-1)
        !       - Kinetic : 0.5 * rho * (vx^2 + vy^2)
        !       - Magnetic: 0.5 * (Bx^2 + By^2)
        halfBL2  = 0.5d0 * (Bx_L*Bx_L + By_L*By_L)
        halfBR2  = 0.5d0 * (Bx_R*Bx_R + By_R*By_R)
        en_L = (P_L-halfBL2)/(gamma-1.d0) + 0.5d0*rho_L*(vx_L*vx_L + vy_L*vy_L) + halfBL2
        en_R = (P_R-halfBR2)/(gamma-1.d0) + 0.5d0*rho_R*(vx_R*vx_R + vy_R*vy_R) + halfBR2

        ! Average states for primitive and conserved variables
        rho_avg  = 0.5d0 * (rho_L + rho_R)
        inv_rho_avg = 1.0d0 / rho_avg
        momx_avg = 0.5d0 * (rho_L*vx_L + rho_R*vx_R)
        momy_avg = 0.5d0 * (rho_L*vy_L + rho_R*vy_R)
        en_avg   = 0.5d0 * (en_L + en_R)
        Bx_avg   = 0.5d0 * (Bx_L + Bx_R)
        By_avg   = 0.5d0 * (By_L + By_R)

        ! Calculate avg pressure using ideal MHD relation
        P_avg = (gamma - 1.d0) * (en_avg &
                - 0.5d0 * (momx_avg*momx_avg + momy_avg*momy_avg) * inv_rho_avg &
                - 0.5d0 * (Bx_avg*Bx_avg     + By_avg*By_avg)) &
                + 0.5d0 * (Bx_avg*Bx_avg     + By_avg*By_avg)

        ! Calculate Lax-Friedrichs (Rusanov) fluxes for each quantity
        flux_Mass   = momx_avg
        flux_Momx   = momx_avg*momx_avg * inv_rho_avg + P_avg - Bx_avg*Bx_avg
        flux_Momy   = momx_avg*momy_avg * inv_rho_avg - Bx_avg*By_avg
        flux_Energy = (en_avg + P_avg)  * momx_avg    * inv_rho_avg &
                    - Bx_avg * (Bx_avg*momx_avg + By_avg*momy_avg)  &
                    * inv_rho_avg
        flux_By     = (By_avg*momx_avg - Bx_avg*momy_avg) * inv_rho_avg

        ! Estimate max local wavespeed
        c_L2 = (gamma * (P_L - halfBL2) + 2.0d0 * halfBL2) / rho_L
        c_R2 = (gamma * (P_R - halfBR2) + 2.0d0 * halfBR2) / rho_R
        C_L = sqrt( 0.5d0 * (c_L2 + abs(c_L2))) + abs(vx_L)
        C_R = sqrt( 0.5d0 * (c_R2 + abs(c_R2))) + abs(vx_R)
        C   = 0.5d0 * max(C_L, C_R) 

        ! Add stabilizing diffusive term (Rusanov/Lax-Friedrichs) to reduce oscillations
        flux_Mass   = flux_Mass   - C * (rho_L      - rho_R)
        flux_Momx   = flux_Momx   - C * (rho_L*vx_L - rho_R*vx_R)
        flux_Momy   = flux_Momy   - C * (rho_L*vy_L - rho_R*vy_R)
        flux_Energy = flux_Energy - C * (en_L       - en_R)
        flux_By     = flux_By     - C * (By_L       - By_R)

    end subroutine compute_fluxes


end module mhd_flux