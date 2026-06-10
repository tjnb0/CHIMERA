module mhd_flux
    use mhd_field_ops
    use mhd_config
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    !---------------------------------------------------------------------------
    ! Purpose: Subroutines for flux computations in MHD solver: apply numerical
    !          fluxes to conserved fields, constrained transport for divergence
    !          control, and flux function evaluation. Supports periodic,
    !          outflow, fixed, and inflow BCs via per-side flags in mhd_config.
    !---------------------------------------------------------------------------


    implicit none
    private
    public update_conserved, constrained_transport, compute_fluxes
    
contains


    subroutine update_conserved(F, flux_F_X, flux_F_Y, dx, dt, nx, ny)
    !
    !   Apply fluxes to conserved field F using BC-aware indexing.
    !   For each cell: F(i,j) -= dtdx * (flux_right(i,j) - flux_left(i,j))
    !                          + dtdx * (flux_top(i,j)   - flux_bottom(i,j))  [signs included]
    !
    !   xhi and yhi boundary fluxes are already encoded in flux_F_X(nx,:) and
    !   flux_F_Y(:,ny) by the outflow ghost state set in reconstruction (step 4).
    !   Only xlo and ylo boundaries need special treatment here.
    !
    !   For outflow at xlo/ylo: no incoming flux from outside the domain.
    !   For periodic at xlo/ylo: incoming flux wraps from the opposite boundary.
    !
    !   Inputs:
    !       - flux_F_X : x-direction flux (nx, ny)
    !       - flux_F_Y : y-direction flux (nx, ny)
    !       - dx       : cell size
    !       - dt       : time step
    !
    !   InOut:
    !       - F : conserved variable (nx, ny)
    !
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: dx, dt
        real(8), intent(in) :: flux_F_X(nx, ny), flux_F_Y(nx, ny)
        real(8), intent(inout) :: F(nx, ny)
        real(8) :: dtdx

        dtdx = dt * dx

        ! --- X-direction update ---

        ! Interior cells (i = 2..nx): incoming flux from left neighbor
        F(2:nx, :) = F(2:nx, :) - dtdx*flux_F_X(2:nx, :) + dtdx*flux_F_X(1:nx-1, :)

        ! xlo boundary (i = 1): incoming left-face flux depends on BC
        select case (bc_xlo)
            case (BC_PERIODIC)
                F(1, :) = F(1, :) - dtdx*flux_F_X(1, :) + dtdx*flux_F_X(nx, :)
            case default  ! outflow/fixed/inflow: no incoming flux from outside
                F(1, :) = F(1, :) - dtdx*flux_F_X(1, :)
        end select

        ! --- Y-direction update ---

        ! Interior cells (j = 2..ny): incoming flux from bottom neighbor
        F(:, 2:ny) = F(:, 2:ny) - dtdx*flux_F_Y(:, 2:ny) + dtdx*flux_F_Y(:, 1:ny-1)

        ! ylo boundary (j = 1): incoming bottom-face flux depends on BC
        select case (bc_ylo)
            case (BC_PERIODIC)
                F(:, 1) = F(:, 1) - dtdx*flux_F_Y(:, 1) + dtdx*flux_F_Y(:, ny)
            case default  ! outflow/fixed/inflow: no incoming flux from outside
                F(:, 1) = F(:, 1) - dtdx*flux_F_Y(:, 1)
        end select

    end subroutine update_conserved


    subroutine constrained_transport(bx, by, flux_By_X, flux_Bx_Y, dx, dt, nx, ny)
    !
    !   Update face-centred B fields via constrained transport (CT).
    !   Ez at each cell corner (i+1/2, j+1/2) is assembled from the four
    !   surrounding cell fluxes. The two cshift calls that gathered the
    !   (i, j+1) and (i+1, j) neighbors are replaced with explicit BC-aware
    !   array slices: interior cells use direct neighbors; boundary cells use
    !   either a periodic wrap or zero-gradient copy of the nearest interior flux.
    !
    !   Inputs:
    !       - flux_By_X : x-dir flux of By (nx, ny)
    !       - flux_Bx_Y : y-dir flux of Bx (nx, ny)
    !       - dx        : cell size
    !       - dt        : time step
    !
    !   InOuts:
    !       - bx : face-centred Bx
    !       - by : face-centred By
    !
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: dx, dt
        real(8), intent(in) :: flux_By_X(nx, ny), flux_Bx_Y(nx, ny)
        real(8), intent(inout) :: bx(nx, ny), by(nx, ny)

        real(8) :: Ez(nx, ny)
        real(8) :: dbx(nx, ny), dby(nx, ny)
        real(8) :: flux_By_X_up(nx, ny)    ! flux_By_X shifted up   (j -> j+1)
        real(8) :: flux_Bx_Y_right(nx, ny) ! flux_Bx_Y shifted right (i -> i+1)

        ! --- Gather (i, j+1) neighbor of flux_By_X ---
        ! Interior: direct slice
        flux_By_X_up(:, 1:ny-1) = flux_By_X(:, 2:ny)
        ! yhi boundary: depends on BC
        select case (bc_yhi)
            case (BC_PERIODIC)
                flux_By_X_up(:, ny) = flux_By_X(:, 1)
            case default  ! zero-gradient: copy nearest interior value
                flux_By_X_up(:, ny) = flux_By_X(:, ny)
        end select

        ! --- Gather (i+1, j) neighbor of flux_Bx_Y ---
        flux_Bx_Y_right(1:nx-1, :) = flux_Bx_Y(2:nx, :)
        select case (bc_xhi)
            case (BC_PERIODIC)
                flux_Bx_Y_right(nx, :) = flux_Bx_Y(1, :)
            case default
                flux_Bx_Y_right(nx, :) = flux_Bx_Y(nx, :)
        end select

        ! --- Ez at cell corners from surrounding face fluxes ---
        Ez = 0.25d0 * (-flux_By_X - flux_By_X_up + flux_Bx_Y + flux_Bx_Y_right)

        ! --- Discrete curl of -Ez updates bx and by ---
        call compute_curl_2d(-Ez, dx, nx, ny, dbx, dby)
        bx = bx + dt * dbx
        by = by + dt * dby

    end subroutine constrained_transport


    subroutine compute_fluxes(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R,   &
                              Bx_L, Bx_R, By_L, By_R, gamma, nx, ny,            &
                              flux_Mass, flux_Momx, flux_Momy, flux_Energy, flux_By)
    !   Routes to rusanov_flux or hlld_flux based on riemann_solver selected
    !   in mhd_config.f90.  .
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
        integer, intent(in)    :: nx, ny
        real(8), intent(inout) :: rho_L(nx, ny), rho_R(nx, ny)
        real(8), intent(in)    :: vx_L(nx, ny), vx_R(nx, ny)
        real(8), intent(in)    :: vy_L(nx, ny), vy_R(nx, ny)
        real(8), intent(inout) :: P_L(nx, ny), P_R(nx, ny)
        real(8), intent(in)    :: Bx_L(nx, ny), Bx_R(nx, ny)
        real(8), intent(in)    :: By_L(nx, ny), By_R(nx, ny)
        real(8), intent(in)    :: gamma
        real(8), intent(out)   :: flux_Mass(nx, ny), flux_Momx(nx, ny), flux_Momy(nx, ny)
        real(8), intent(out)   :: flux_Energy(nx, ny), flux_By(nx, ny)

        select case (riemann_solver)
        case (RIEMANN_RUSANOV)
            call rusanov_flux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R, &
                              Bx_L, Bx_R, By_L, By_R, gamma, nx, ny,           &
                              flux_Mass, flux_Momx, flux_Momy, flux_Energy, flux_By)
        case (RIEMANN_HLLE)
            call hlle_flux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R,    &
                           Bx_L, Bx_R, By_L, By_R, gamma, nx, ny,              &
                           flux_Mass, flux_Momx, flux_Momy, flux_Energy, flux_By)
        case (RIEMANN_HLLD)
            call hlld_flux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R,    &
                           Bx_L, Bx_R, By_L, By_R, gamma, nx, ny,              &
                           flux_Mass, flux_Momx, flux_Momy, flux_Energy, flux_By)
        end select

    end subroutine compute_fluxes


    ! =========================================================================
    ! Private Riemann solvers
    ! =========================================================================
    subroutine rusanov_flux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R,   &
                            Bx_L, Bx_R, By_L, By_R, gamma, nx, ny,             &
                            flux_Mass, flux_Momx, flux_Momy, flux_Energy, flux_By)
    !   Local Lax-Friedrichs (Rusanov) flux.
    !
    !   Method:
    !       1. Compute energies of left and right states.
    !       2. Compute averaged states of primitive and conserved quantities.
    !       3. Evaluate fluxes from the star states.
    !       4. Estimate local maximum wavespeed using fast magnetosonic speed.
    !       5. Apply a stabilizing diffusive term proportional to the wavespeed.
    !
    !   Notes:
    !       - Assumes ideal MHD equations in 2D.
    !       - Magnetic field is split into Bx (normal) and By (tangential) components.
    !         main.f90 passes rotated arguments for y-direction faces.
    !
        integer, intent(in)    :: nx, ny
        real(8), intent(inout) :: rho_L(nx, ny), rho_R(nx, ny)
        real(8), intent(in)    :: vx_L(nx, ny), vx_R(nx, ny)
        real(8), intent(in)    :: vy_L(nx, ny), vy_R(nx, ny)
        real(8), intent(inout) :: P_L(nx, ny), P_R(nx, ny)
        real(8), intent(in)    :: Bx_L(nx, ny), Bx_R(nx, ny)
        real(8), intent(in)    :: By_L(nx, ny), By_R(nx, ny)
        real(8), intent(in)    :: gamma
        real(8), intent(out)   :: flux_Mass(nx, ny), flux_Momx(nx, ny), flux_Momy(nx, ny)
        real(8), intent(out)   :: flux_Energy(nx, ny), flux_By(nx, ny)

        ! Local arrays
        real(8) :: en_L(nx, ny), en_R(nx, ny), halfBL2(nx, ny), halfBR2(nx, ny)
        real(8) :: rho_avg(nx, ny), inv_rho_avg(nx, ny)
        real(8) :: momx_avg(nx, ny), momy_avg(nx, ny), en_avg(nx, ny)
        real(8) :: Bx_avg(nx, ny), By_avg(nx, ny), P_avg(nx, ny)
        real(8) :: c_L2(nx, ny), c_R2(nx, ny), C_L(nx, ny), C_R(nx, ny), C(nx, ny)
        real(8) :: p_th_L(nx, ny), p_th_R(nx, ny)

        ! Compute total energy for left and right states
        !   - en = internal + kinetic + magnetic energy
        !       - Internal: (P - 0.5 * B^2) / (gamma-1)
        !       - Kinetic : 0.5 * rho * (vx^2 + vy^2)
        !       - Magnetic: 0.5 * (Bx^2 + By^2)
        halfBL2  = 0.5d0 * (Bx_L*Bx_L + By_L*By_L)
        halfBR2  = 0.5d0 * (Bx_R*Bx_R + By_R*By_R)

        ! Compute floored thermal pressures locally to prevent negative
        ! thermal pressure during energy or underestimating wave speed
        p_th_L = max(P_L - halfBL2, P_floor)
        p_th_R = max(P_R - halfBR2, P_floor)

        en_L = p_th_L/(gamma-1.d0) + 0.5d0*rho_L*(vx_L*vx_L + vy_L*vy_L) + halfBL2
        en_R = p_th_R/(gamma-1.d0) + 0.5d0*rho_R*(vx_R*vx_R + vy_R*vy_R) + halfBR2

        ! Average states for primitive and conserved variables
        rho_avg     = 0.5d0 * (rho_L + rho_R)
        inv_rho_avg = 1.0d0 / rho_avg
        momx_avg    = 0.5d0 * (rho_L*vx_L + rho_R*vx_R)
        momy_avg    = 0.5d0 * (rho_L*vy_L + rho_R*vy_R)
        en_avg      = 0.5d0 * (en_L + en_R)
        Bx_avg      = 0.5d0 * (Bx_L + Bx_R)
        By_avg      = 0.5d0 * (By_L + By_R)

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

        ! Estimate max local wavespeed using floored thermal pressures
        ! so a near-zero or negative p_thermal cannot underestimate c_f
        c_L2 = (gamma * p_th_L + 2.0d0 * halfBL2) / rho_L
        c_R2 = (gamma * p_th_R + 2.0d0 * halfBR2) / rho_R
        C_L = sqrt( 0.5d0 * (c_L2 + abs(c_L2))) + abs(vx_L)
        C_R = sqrt( 0.5d0 * (c_R2 + abs(c_R2))) + abs(vx_R)
        C   = 0.5d0 * max(C_L, C_R)

        ! Add stabilizing diffusive term (Rusanov/Lax-Friedrichs) to reduce oscillations
        flux_Mass   = flux_Mass   - C * (rho_L      - rho_R)
        flux_Momx   = flux_Momx   - C * (rho_L*vx_L - rho_R*vx_R)
        flux_Momy   = flux_Momy   - C * (rho_L*vy_L - rho_R*vy_R)
        flux_Energy = flux_Energy - C * (en_L       - en_R)
        flux_By     = flux_By     - C * (By_L       - By_R)

    end subroutine rusanov_flux


    subroutine hlle_flux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R,   &
                     Bx_L, Bx_R, By_L, By_R, gamma, nx, ny,             &
                     flux_Mass, flux_Momx, flux_Momy, flux_Energy, flux_By)
    !   HLLE (Harten-Lax-van Leer-Einfeldt) flux for ideal 2D MHD.
    !
    !   Method:
    !       1. Compute floored thermal pressures and total energies for L/R states.
    !       2. Estimate fast magnetosonic speed c_f for each state using the same
    !          upper bound as Rusanov: c_f^2 = (gamma*p_th + 2*halfB2) / rho.
    !       3. Davis-type signal speed estimates:
    !              S_L = min(vx_L - c_fL, vx_R - c_fR)
    !              S_R = max(vx_L + c_fL, vx_R + c_fR)
    !       4. Evaluate physical fluxes F_L, F_R from each state directly.
    !       5. Assemble HLL intercell flux:
    !              F = (S_R*F_L - S_L*F_R + S_L*S_R*(U_R - U_L)) / (S_R - S_L)
    !          with supersonic limiting:
    !              F = F_L  if S_L >= 0  (all waves right-going)
    !              F = F_R  if S_R <= 0  (all waves left-going)
    !
    !   Notes:
    !       - P_L, P_R are TOTAL pressures (thermal + magnetic), consistent with
    !         get_primitive and rusanov_flux.
    !       - Using the same c_f upper bound as Rusanov keeps the two solvers
    !         directly comparable and avoids disagreements in the CFL estimate.
    !       - main.f90 passes rotated arguments for y-direction faces.
    !
        integer, intent(in)    :: nx, ny
        real(8), intent(inout) :: rho_L(nx, ny), rho_R(nx, ny)
        real(8), intent(in)    :: vx_L(nx, ny),  vx_R(nx, ny)
        real(8), intent(in)    :: vy_L(nx, ny),  vy_R(nx, ny)
        real(8), intent(inout) :: P_L(nx, ny),   P_R(nx, ny)
        real(8), intent(in)    :: Bx_L(nx, ny),  Bx_R(nx, ny)
        real(8), intent(in)    :: By_L(nx, ny),  By_R(nx, ny)
        real(8), intent(in)    :: gamma
        real(8), intent(out)   :: flux_Mass(nx, ny), flux_Momx(nx, ny), flux_Momy(nx, ny)
        real(8), intent(out)   :: flux_Energy(nx, ny), flux_By(nx, ny)

        ! Local arrays
        real(8) :: halfBL2(nx, ny), halfBR2(nx, ny)
        real(8) :: p_th_L(nx, ny), p_th_R(nx, ny)
        real(8) :: en_L(nx, ny), en_R(nx, ny)
        real(8) :: c_fL(nx, ny), c_fR(nx, ny)
        real(8) :: S_L(nx, ny), S_R(nx, ny), inv_Sdiff(nx, ny)
        ! Left-state physical fluxes
        real(8) :: FL_Mass(nx, ny), FL_Momx(nx, ny), FL_Momy(nx, ny)
        real(8) :: FL_Energy(nx, ny), FL_By(nx, ny)
        ! Right-state physical fluxes
        real(8) :: FR_Mass(nx, ny), FR_Momx(nx, ny), FR_Momy(nx, ny)
        real(8) :: FR_Energy(nx, ny), FR_By(nx, ny)

        ! --- Step 1: Magnetic energy and floored thermal pressures ---
        halfBL2 = 0.5d0 * (Bx_L*Bx_L + By_L*By_L)
        halfBR2 = 0.5d0 * (Bx_R*Bx_R + By_R*By_R)
        p_th_L  = max(P_L - halfBL2, P_floor)
        p_th_R  = max(P_R - halfBR2, P_floor)

        ! --- Step 2: Total energies ---
        en_L = p_th_L/(gamma - 1.d0) + 0.5d0*rho_L*(vx_L*vx_L + vy_L*vy_L) + halfBL2
        en_R = p_th_R/(gamma - 1.d0) + 0.5d0*rho_R*(vx_R*vx_R + vy_R*vy_R) + halfBR2

        ! --- Step 3: Fast magnetosonic speeds ---
        ! c_f^2 = (gamma*p_th + 2*halfB2) / rho; safe sqrt via 0.5*(x+|x|) = max(x,0)
        c_fL = sqrt(0.5d0 * ((gamma*p_th_L + 2.d0*halfBL2) / rho_L &
                            + abs((gamma*p_th_L + 2.d0*halfBL2) / rho_L)))
        c_fR = sqrt(0.5d0 * ((gamma*p_th_R + 2.d0*halfBR2) / rho_R &
                            + abs((gamma*p_th_R + 2.d0*halfBR2) / rho_R)))

        ! --- Step 4: Davis-type signal speed estimates ---
        S_L = min(vx_L - c_fL, vx_R - c_fR)
        S_R = max(vx_L + c_fL, vx_R + c_fR)

        ! --- Step 5: Left-state physical fluxes ---
        FL_Mass   = rho_L * vx_L
        FL_Momx   = rho_L*vx_L*vx_L + P_L - Bx_L*Bx_L
        FL_Momy   = rho_L*vx_L*vy_L - Bx_L*By_L
        FL_Energy = (en_L + P_L)*vx_L - Bx_L*(Bx_L*vx_L + By_L*vy_L)
        FL_By     = By_L*vx_L - Bx_L*vy_L

        ! --- Step 6: Right-state physical fluxes ---
        FR_Mass   = rho_R * vx_R
        FR_Momx   = rho_R*vx_R*vx_R + P_R - Bx_R*Bx_R
        FR_Momy   = rho_R*vx_R*vy_R - Bx_R*By_R
        FR_Energy = (en_R + P_R)*vx_R - Bx_R*(Bx_R*vx_R + By_R*vy_R)
        FR_By     = By_R*vx_R - Bx_R*vy_R

        ! --- Step 7: HLL intercell flux ---
        ! F = (S_R*F_L - S_L*F_R + S_L*S_R*(U_R - U_L)) / (S_R - S_L)
        ! S_R - S_L >= 0 by construction; floor guards the degenerate equal-speed case.
        inv_Sdiff = 1.d0 / max(S_R - S_L, 1.d-12)
        flux_Mass   = (S_R*FR_Mass   - S_L*FL_Mass   + S_L*S_R*(rho_L      - rho_R)     ) * inv_Sdiff
        flux_Momx   = (S_R*FR_Momx   - S_L*FL_Momx   + S_L*S_R*(rho_L*vx_L - rho_R*vx_R)) * inv_Sdiff
        flux_Momy   = (S_R*FR_Momy   - S_L*FL_Momy   + S_L*S_R*(rho_L*vy_L - rho_R*vy_R)) * inv_Sdiff
        flux_Energy = (S_R*FR_Energy - S_L*FL_Energy + S_L*S_R*(en_L       - en_R)      ) * inv_Sdiff
        flux_By     = (S_R*FR_By     - S_L*FL_By     + S_L*S_R*(By_L       - By_R)      ) * inv_Sdiff

        ! --- Step 8: Supersonic limiting ---
        ! Overwrite with exact upwind flux when the Riemann fan doesn't straddle
        ! the interface. The HLL formula is not exact in these limits due to the
        ! inv_Sdiff guard, so explicit where-blocks are required.
        where (S_L >= 0.d0)
            flux_Mass = FR_Mass;  flux_Momx = FR_Momx
            flux_Momy = FR_Momy;  flux_Energy = FR_Energy;  flux_By = FR_By
        end where
        where (S_R <= 0.d0)
            flux_Mass = FL_Mass;  flux_Momx = FL_Momx
            flux_Momy = FL_Momy;  flux_Energy = FL_Energy;  flux_By = FL_By
        end where

    end subroutine hlle_flux


    !subroutine hlle_flux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R,   &
    !                 Bx_L, Bx_R, By_L, By_R, gamma, nx, ny,             &
    !                 flux_Mass, flux_Momx, flux_Momy, flux_Energy, flux_By)
    !end subroutine hlle_flux

end module mhd_flux
