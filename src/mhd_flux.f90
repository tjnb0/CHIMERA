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
            ! Rusanov
            call rusanov_flux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R, &
                              Bx_L, Bx_R, By_L, By_R, gamma, nx, ny,           &
                              flux_Mass, flux_Momx, flux_Momy, flux_Energy, flux_By)
        case (RIEMANN_HLLE)
            ! HLLE
            call hlle_flux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R,    &
                           Bx_L, Bx_R, By_L, By_R, gamma, nx, ny,              &
                           flux_Mass, flux_Momx, flux_Momy, flux_Energy, flux_By)
        case (RIEMANN_HLLD)
            ! HLLD
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

        ! Add stabilizing diffusive term (Rusanov/Lax-Friedrichs) to reduce oscillations.
        flux_Mass   = flux_Mass   - C * (rho_R      - rho_L)
        flux_Momx   = flux_Momx   - C * (rho_R*vx_R - rho_L*vx_L)
        flux_Momy   = flux_Momy   - C * (rho_R*vy_R - rho_L*vy_L)
        flux_Energy = flux_Energy - C * (en_R        - en_L)
        flux_By     = flux_By     - C * (By_R        - By_L)

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
        ! Stability
        real(8) :: cf_max_hlle = 1.0d2

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
        c_fL = min(c_fL, cf_max_hlle)
        c_fR = min(c_fR, cf_max_hlle)

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
        flux_Mass   = (S_R*FL_Mass   - S_L*FR_Mass   + S_L*S_R*(rho_R      - rho_L)     ) * inv_Sdiff
        flux_Momx   = (S_R*FL_Momx   - S_L*FR_Momx   + S_L*S_R*(rho_R*vx_R - rho_L*vx_L)) * inv_Sdiff
        flux_Momy   = (S_R*FL_Momy   - S_L*FR_Momy   + S_L*S_R*(rho_R*vy_R - rho_L*vy_L)) * inv_Sdiff
        flux_Energy = (S_R*FL_Energy - S_L*FR_Energy  + S_L*S_R*(en_R       - en_L)      ) * inv_Sdiff
        flux_By     = (S_R*FL_By     - S_L*FR_By      + S_L*S_R*(By_R       - By_L)      ) * inv_Sdiff

        ! --- Step 8: Supersonic limiting ---
        ! Overwrite with exact upwind flux when the Riemann fan doesn't straddle
        ! the interface. The HLL formula is not exact in these limits due to the
        ! inv_Sdiff guard, so explicit where-blocks are required.
        where (S_L >= 0.d0)                ! All waves right-going: use left state
            flux_Mass   = FL_Mass;   flux_Momx   = FL_Momx
            flux_Momy   = FL_Momy;   flux_Energy = FL_Energy
            flux_By     = FL_By
        end where
        where (S_R <= 0.d0)                ! All waves left-going: use right state
            flux_Mass   = FR_Mass;   flux_Momx   = FR_Momx
            flux_Momy   = FR_Momy;   flux_Energy = FR_Energy
            flux_By     = FR_By
        end where

    end subroutine hlle_flux


    subroutine hlld_flux(rho_L, rho_R, vx_L, vx_R, vy_L, vy_R, P_L, P_R,    &
                         Bx_L, Bx_R, By_L, By_R, gamma, nx, ny,              &
                         flux_Mass, flux_Momx, flux_Momy, flux_Energy, flux_By)
    !   HLLD Riemann solver -- Miyoshi & Kusano (2005), J. Comput. Phys. 208, 315-344.
    !
    !   Resolves seven MHD waves (2 fast, 2 Alfven, 2 slow, 1 entropy) for less
    !   numerical diffusion than Rusanov, especially in current sheets.
    !
    !   Arguments:
    !       vx / Bx = normal  (perpendicular to face)
    !       vy / By = tangential
    !   main.f90 passes rotated arguments for y-direction faces.
    !
    !   Degenerate cases handled:
    !       |DL| or |DR| ~ 0   -> tangential quantities unchanged from L/R state
    !       |Bx| ~ 0           -> Alfven speeds collapse to SM; double-star = outer star
    !       SL >= 0 or SR <= 0 -> supersonic: pure upwind
    !
    !   Reference equations are tagged (M&K Eq. N) in comments below.
    !
        integer, intent(in)    :: nx, ny
        real(8), intent(inout) :: rho_L(nx, ny), rho_R(nx, ny)
        real(8), intent(in)    :: vx_L(nx, ny),  vx_R(nx, ny)
        real(8), intent(in)    :: vy_L(nx, ny),  vy_R(nx, ny)
        real(8), intent(inout) :: P_L(nx, ny),   P_R(nx, ny)
        real(8), intent(in)    :: Bx_L(nx, ny),  Bx_R(nx, ny)
        real(8), intent(in)    :: By_L(nx, ny),  By_R(nx, ny)
        real(8), intent(in)    :: gamma
        real(8), intent(out)   :: flux_Mass(nx, ny), flux_Momx(nx, ny)
        real(8), intent(out)   :: flux_Momy(nx, ny), flux_Energy(nx, ny)
        real(8), intent(out)   :: flux_By(nx, ny)

        integer :: i, j

        ! --- scalar temporaries (one interface at a time) ---
        real(8) :: rhoL, rhoR, vxL, vxR, vyL, vyR, pTL, pTR, BxL, BxR, ByL, ByR
        real(8) :: halfBL2, halfBR2, p_thL, p_thR, enL, enR
        real(8) :: vdotBL, vdotBR

        ! Physical (non-diffusive) fluxes from L and R states
        real(8) :: fM_L, fMx_L, fMy_L, fE_L, fB_L
        real(8) :: fM_R, fMx_R, fMy_R, fE_R, fB_R

        ! Fast-wave and outer wave speeds
        real(8) :: csL2, csR2, b2L, b2R, caL, caR, tmpL, tmpR, cfL, cfR
        real(8) :: SL, SR

        ! Contact speed and star pressure
        real(8) :: SM, pT_star, denMD

        ! HLL fallback (Prong 2: SM outside [SL,SR])
        real(8) :: inv_SRSL

        ! Minimum margin (as a fraction of fan width SR-SL) that the contact speed
        ! SM must keep from each fan edge. If SM is closer than EPS_FAN*(SR-SL) to
        ! SL or SR, the star-state denominators (S-SM) are ill-conditioned and HLL
        ! is used instead. Matches the eps = 1e-6 degeneracy factor of
        ! Matsumoto, Miyoshi & Takasao (2019), Sec. 2.6.
        real(8), parameter :: EPS_FAN = 1.0d-6

        ! Left outer star state
        real(8) :: rho_Ls, denL, vy_Ls, By_Ls, vdotBLs, en_Ls

        ! Right outer star state
        real(8) :: rho_Rs, denR, vy_Rs, By_Rs, vdotBRs, en_Rs

        ! Alfven speeds and double-star (inner) states
        real(8) :: SqrtRhoLs, SqrtRhoRs, Bx_avg, signBx
        real(8) :: SLs, SRs, denom_ss
        real(8) :: vy_ss, By_ss, vdotBss, en_Lss, en_Rss

        ! Catch-all robustness guard: HLLD flux temporaries, bounds, and HLL fallback
        real(8) :: fmass, fmomx, fmomy, fener, fbyf
        real(8) :: Smax, bnd_mass, bnd_momx, bnd_momy, bnd_enr, bnd_by
        logical :: use_hll

        ! Acceptance factor for the catch-all bound. A consistent Godunov flux is
        ! bounded by max|F_L,F_R| + Smax*(|U_L|+|U_R|); FAC_BND gives generous slack
        ! so only genuinely degenerate cells fall back to HLL.
        real(8), parameter :: FAC_BND = 8.0d0
        real(8), parameter :: cf_max_hlld = 50.0d0

        do j = 1, ny
            do i = 1, nx

                ! --- extract left / right primitive states ---
                rhoL = rho_L(i,j); rhoR = rho_R(i,j)
                vxL  = vx_L(i,j);  vxR  = vx_R(i,j)
                vyL  = vy_L(i,j);  vyR  = vy_R(i,j)
                pTL  = P_L(i,j);   pTR  = P_R(i,j)    ! total pressure
                BxL  = Bx_L(i,j);  BxR  = Bx_R(i,j)
                ByL  = By_L(i,j);  ByR  = By_R(i,j)

                ! floored thermal pressures and total energies
                halfBL2 = 0.5d0*(BxL**2 + ByL**2)
                halfBR2 = 0.5d0*(BxR**2 + ByR**2)
                p_thL   = max(pTL - halfBL2, P_floor)
                p_thR   = max(pTR - halfBR2, P_floor)
                enL = p_thL/(gamma-1.0d0) + 0.5d0*rhoL*(vxL**2 + vyL**2) + halfBL2
                enR = p_thR/(gamma-1.0d0) + 0.5d0*rhoR*(vxR**2 + vyR**2) + halfBR2

                ! physical MHD fluxes in the normal (x) direction
                vdotBL = vxL*BxL + vyL*ByL
                fM_L   = rhoL*vxL
                fMx_L  = rhoL*vxL**2  + pTL - BxL**2
                fMy_L  = rhoL*vxL*vyL - BxL*ByL
                fE_L   = (enL + pTL)*vxL - BxL*vdotBL
                fB_L   = ByL*vxL - BxL*vyL

                vdotBR = vxR*BxR + vyR*ByR
                fM_R   = rhoR*vxR
                fMx_R  = rhoR*vxR**2  + pTR - BxR**2
                fMy_R  = rhoR*vxR*vyR - BxR*ByR
                fE_R   = (enR + pTR)*vxR - BxR*vdotBR
                fB_R   = ByR*vxR - BxR*vyR

                ! --- fast magnetosonic speeds ---
                caL  = BxL**2 / rhoL              ! normal Alfven speed^2
                caR  = BxR**2 / rhoR
                csL2 = gamma*p_thL / rhoL         ! sound speed^2
                csR2 = gamma*p_thR / rhoR
                b2L  = (BxL**2 + ByL**2) / rhoL   ! total Alfven speed^2
                b2R  = (BxR**2 + ByR**2) / rhoR

                ! cf = sqrt(0.5*(cs^2 + b^2 + sqrt((cs^2+b^2)^2 - 4*cs^2*ca^2)))
                tmpL = csL2 + b2L
                tmpR = csR2 + b2R
                cfL  = sqrt(0.5d0*(tmpL + sqrt(max(tmpL**2 - 4.0d0*csL2*caL, 0.0d0))))
                cfR  = sqrt(0.5d0*(tmpR + sqrt(max(tmpR**2 - 4.0d0*csR2*caR, 0.0d0))))
                cfL = min(cfL, cf_max_hlld)
                cfR = min(cfR, cf_max_hlld)

                ! outer wave speeds: Davis (1988) estimates
                SL = min(vxL - cfL, vxR - cfR)
                SR = max(vxL + cfL, vxR + cfR)

                ! Safety catch 1: 
                ! ---------------
                ! Unphysical wave speed from near-vacuum reconstruction or bad state.
                ! Use upwind from the slower side.
                if (cfL >= cf_max_hlld .or. cfR >= cf_max_hlld) then
                    if (cfL <= cfR) then
                        flux_Mass(i,j)   = fM_L;   flux_Momx(i,j)   = fMx_L
                        flux_Momy(i,j)   = fMy_L;  flux_Energy(i,j) = fE_L
                        flux_By(i,j)     = fB_L
                    else
                        flux_Mass(i,j)   = fM_R;   flux_Momx(i,j)   = fMx_R
                        flux_Momy(i,j)   = fMy_R;  flux_Energy(i,j) = fE_R
                        flux_By(i,j)     = fB_R
                    end if
                    cycle
                end if

                ! --- supersonic cases: upwind ---
                if (SL >= 0.0d0) then
                    flux_Mass(i,j)   = fM_L;   flux_Momx(i,j)   = fMx_L
                    flux_Momy(i,j)   = fMy_L;  flux_Energy(i,j) = fE_L
                    flux_By(i,j)     = fB_L
                    cycle
                end if
                if (SR <= 0.0d0) then
                    flux_Mass(i,j)   = fM_R;   flux_Momx(i,j)   = fMx_R
                    flux_Momy(i,j)   = fMy_R;  flux_Energy(i,j) = fE_R
                    flux_By(i,j)     = fB_R
                    cycle
                end if

                ! --- contact speed SM and star total pressure pT* ---
                ! (M&K Eq. 38; pT* from Rankine-Hugoniot across outer waves)
                denMD = rhoR*(SR - vxR) - rhoL*(SL - vxL)
                if (abs(denMD) < 1.0d-15) &
                    denMD = sign(1.0d-15, denMD)

                SM      = (rhoR*vxR*(SR-vxR) - rhoL*vxL*(SL-vxL) - (pTR - pTL)) / denMD
                pT_star = pTL + rhoL*(SL - vxL)*(SM - vxL)


                ! Safety catch 2: 
                ! ---------------
                ! SM not safely inside the fan -> use HLL
                !
                ! Require SM to be at least eps_fan = EPS_FAN*(SR-SL) away from each
                ! edge. This bounds every (S-SM) denominator >= eps_fan > 0, keeping
                ! all star quantities finite. An extra HLL diffusion should be small.
                ! cf. Matsumoto, Miyoshi & Takasao (2019), Sec. 2.6 (the eps(SR-SL)
                ! denominator-degeneracy guard).
                if (SM <= SL + EPS_FAN*(SR-SL) .or. SM >= SR - EPS_FAN*(SR-SL)) then
                    inv_SRSL         = 1.0d0 / (SR - SL)
                    flux_Mass(i,j)   = (SR*fM_L  - SL*fM_R  + SL*SR*(rhoR     - rhoL    )) * inv_SRSL
                    flux_Momx(i,j)   = (SR*fMx_L - SL*fMx_R + SL*SR*(rhoR*vxR - rhoL*vxL)) * inv_SRSL
                    flux_Momy(i,j)   = (SR*fMy_L - SL*fMy_R + SL*SR*(rhoR*vyR - rhoL*vyL)) * inv_SRSL
                    flux_Energy(i,j) = (SR*fE_L  - SL*fE_R  + SL*SR*(enR      - enL     )) * inv_SRSL
                    flux_By(i,j)     = (SR*fB_L  - SL*fB_R  + SL*SR*(ByR      - ByL     )) * inv_SRSL
                    cycle
                end if

                ! Safety catch 3: 
                ! ---------------
                ! pT_star slightly negative from Davis wave-speed under-estimate.
                ! Clamp to P_floor. 
                pT_star = max(pT_star, P_floor)

                ! --- left outer star state (M&K Eqs. 22-23, 28-29, 31) ---
                rho_Ls = rhoL*(SL - vxL) / (SL - SM)
                denL   = rhoL*(SL - vxL)*(SL - SM) - BxL**2

                ! denL -> 0 is the S_L = S_L* merge (fast wave = Alfven wave). denL is
                ! the difference of two O(1) terms [rhoL(SL-vxL)(SL-SM)] and [BxL**2];
                ! when they cancel, the vy*/By* corrections below blow up. Compare
                ! |denL| to the magnitude of BOTH terms (not just one) so partial
                ! cancellation is caught and the tangential fields fall back to the
                ! L state (correct in the merged-wave limit).
                if (abs(denL) > 1.0d-3 * max(abs(rhoL*(SL-vxL)*(SL-SM)), BxL**2, 1.0d-30)) then
                    vy_Ls = vyL - BxL*ByL*(SM  - vxL) / denL              ! Eq. 28
                    By_Ls = ByL*(rhoL*(SL-vxL)**2 - BxL**2) / denL        ! Eq. 29
                else
                    vy_Ls = vyL;   By_Ls = ByL    ! degenerate: fast wave = Alfven wave
                end if

                vdotBLs = SM*BxL + vy_Ls*By_Ls
                en_Ls   = (enL*(SL-vxL) - pTL*vxL + pT_star*SM &
                           + BxL*(vdotBL - vdotBLs)) / (SL - SM)          ! Eq. 31

                ! --- right outer star state ---
                rho_Rs = rhoR*(SR - vxR) / (SR - SM)
                denR   = rhoR*(SR - vxR)*(SR - SM) - BxR**2

                if (abs(denR) > 1.0d-3 * max(abs(rhoR*(SR-vxR)*(SR-SM)), BxR**2, 1.0d-30)) then
                    vy_Rs = vyR - BxR*ByR*(SM  - vxR) / denR
                    By_Rs = ByR*(rhoR*(SR-vxR)**2 - BxR**2) / denR
                else
                    vy_Rs = vyR;   By_Rs = ByR
                end if

                vdotBRs = SM*BxR + vy_Rs*By_Rs
                en_Rs   = (enR*(SR-vxR) - pTR*vxR + pT_star*SM &
                           + BxR*(vdotBR - vdotBRs)) / (SR - SM)

                ! --- Alfven speeds and double-star (inner) states ---
                SqrtRhoLs = sqrt(max(rho_Ls, rho_floor))
                SqrtRhoRs = sqrt(max(rho_Rs, rho_floor))
                Bx_avg    = 0.5d0*(BxL + BxR)    ! continuous across contact (div B = 0)
                SLs       = SM - abs(Bx_avg) / SqrtRhoLs
                SRs       = SM + abs(Bx_avg) / SqrtRhoRs

                if (abs(Bx_avg) > 1.0d-10) then
                    ! General: distinct Alfven waves on each side
                    signBx   = sign(1.0d0, Bx_avg)
                    denom_ss = SqrtRhoLs + SqrtRhoRs

                    vy_ss = (SqrtRhoLs*vy_Ls + SqrtRhoRs*vy_Rs &
                             + signBx*(By_Rs - By_Ls)) / denom_ss          ! Eq. 47
                    By_ss = (SqrtRhoLs*By_Rs + SqrtRhoRs*By_Ls &
                             + signBx*SqrtRhoLs*SqrtRhoRs*(vy_Rs - vy_Ls)) / denom_ss  ! Eq. 49

                    vdotBss = SM*Bx_avg + vy_ss*By_ss
                    en_Lss  = en_Ls - SqrtRhoLs*signBx*(vdotBLs - vdotBss)! Eq. 51
                    en_Rss  = en_Rs + SqrtRhoRs*signBx*(vdotBRs - vdotBss)! Eq. 52
                else
                    ! Degenerate: Bx ~ 0; Alfven speeds collapse onto SM
                    ! Double-star states = average of outer star tangential values
                    vy_ss  = 0.5d0*(vy_Ls + vy_Rs)
                    By_ss  = 0.5d0*(By_Ls + By_Rs)
                    en_Lss = en_Ls
                    en_Rss = en_Rs
                    SLs    = SM
                    SRs    = SM
                end if

                ! --- flux region selection ---
                !
                !  SL     SL*   SM   SR*    SR
                !  |---II--|--III-|--IV--|---V---|
                !  F_L*  F_L** F_R**  F_R*
                !
                ! Flux from any star state U* : F*  = F  + S_outer*(U* - U)
                ! Flux from double-star    U**: F** = F* + S_Alfven*(U** - U*)
                !
                ! Density and normal momentum are unchanged across the Alfven wave
                ! (rho** = rho*, vx** = SM), so those two rows simplify.
                !
                if (SM >= 0.0d0) then
                    if (SLs >= 0.0d0) then
                        ! Region II: left outer star
                        fmass = fM_L  + SL*(rho_Ls       - rhoL)
                        fmomx = fMx_L + SL*(rho_Ls*SM    - rhoL*vxL)
                        fmomy = fMy_L + SL*(rho_Ls*vy_Ls - rhoL*vyL)
                        fener = fE_L  + SL*(en_Ls         - enL)
                        fbyf  = fB_L  + SL*(By_Ls         - ByL)
                    else
                        ! Region III: left double-star (inner Alfven state)
                        fmass = fM_L  + SL*(rho_Ls       - rhoL)
                        fmomx = fMx_L + SL*(rho_Ls*SM    - rhoL*vxL)
                        fmomy = fMy_L + SL*(rho_Ls*vy_Ls - rhoL*vyL) &
                                                  + SLs*rho_Ls*(vy_ss - vy_Ls)
                        fener = fE_L  + SL*(en_Ls         - enL) &
                                                  + SLs*(en_Lss       - en_Ls)
                        fbyf  = fB_L  + SL*(By_Ls         - ByL) &
                                                  + SLs*(By_ss        - By_Ls)
                    end if
                else
                    if (SRs <= 0.0d0) then
                        ! Region V: right outer star
                        fmass = fM_R  + SR*(rho_Rs       - rhoR)
                        fmomx = fMx_R + SR*(rho_Rs*SM    - rhoR*vxR)
                        fmomy = fMy_R + SR*(rho_Rs*vy_Rs - rhoR*vyR)
                        fener = fE_R  + SR*(en_Rs         - enR)
                        fbyf  = fB_R  + SR*(By_Rs         - ByR)
                    else
                        ! Region IV: right double-star (inner Alfven state)
                        fmass = fM_R  + SR*(rho_Rs       - rhoR)
                        fmomx = fMx_R + SR*(rho_Rs*SM    - rhoR*vxR)
                        fmomy = fMy_R + SR*(rho_Rs*vy_Rs - rhoR*vyR) &
                                                  + SRs*rho_Rs*(vy_ss - vy_Rs)
                        fener = fE_R  + SR*(en_Rs         - enR) &
                                                  + SRs*(en_Rss       - en_Rs)
                        fbyf  = fB_R  + SR*(By_Rs         - ByR) &
                                                  + SRs*(By_ss        - By_Rs)
                    end if
                end if

                ! ---------------------------------------------------------------
                ! Catch-all robustness guard.
                !
                ! The five-wave HLLD state involves several divisions (Eqs. 20,
                ! 24-25, 28) whose denominators can become small in degenerate
                ! configurations -- (S-SM)->0, denL/denR->0 (fast=Alfven merge),
                ! sqrt(rho*)->0, etc. Even with the targeted guards above, a
                ! pathological reconstructed state can drive an HLLD flux component
                ! far outside the physically admissible range, injecting a large
                ! flux that seeds a runaway cascade (observed on Orszag-Tang).
                !
                ! Any consistent Godunov flux is bounded by the input physical
                ! fluxes plus the wave-speed times the state jump. We accept the
                ! HLLD flux only if every component is finite and within FAC times
                ! that bound; otherwise we fall back to the (always-bounded) HLL
                ! flux for this interface. The bound is generous (FAC=8), so this
                ! triggers only on genuinely degenerate cells and leaves the HLLD
                ! solution untouched everywhere else.
                Smax = max(abs(SL), abs(SR))
                bnd_mass = max(abs(fM_L),  abs(fM_R))  + Smax*(rhoL      + rhoR)
                bnd_momx = max(abs(fMx_L), abs(fMx_R)) + Smax*(abs(rhoL*vxL) + abs(rhoR*vxR))
                bnd_momy = max(abs(fMy_L), abs(fMy_R)) + Smax*(abs(rhoL*vyL) + abs(rhoR*vyR))
                bnd_enr  = max(abs(fE_L),  abs(fE_R))  + Smax*(abs(enL)   + abs(enR))
                bnd_by   = max(abs(fB_L),  abs(fB_R))  + Smax*(abs(ByL)   + abs(ByR))

                use_hll = .false.
                if (.not. ieee_is_finite(fmass) .or. abs(fmass) > FAC_BND*bnd_mass) use_hll = .true.
                if (.not. ieee_is_finite(fmomx) .or. abs(fmomx) > FAC_BND*bnd_momx) use_hll = .true.
                if (.not. ieee_is_finite(fmomy) .or. abs(fmomy) > FAC_BND*bnd_momy) use_hll = .true.
                if (.not. ieee_is_finite(fener) .or. abs(fener) > FAC_BND*bnd_enr)  use_hll = .true.
                if (.not. ieee_is_finite(fbyf)  .or. abs(fbyf)  > FAC_BND*bnd_by)   use_hll = .true.

                if (use_hll) then
                    inv_SRSL = 1.0d0 / (SR - SL)
                    fmass = (SR*fM_L  - SL*fM_R  + SL*SR*(rhoR     - rhoL    )) * inv_SRSL
                    fmomx = (SR*fMx_L - SL*fMx_R + SL*SR*(rhoR*vxR - rhoL*vxL)) * inv_SRSL
                    fmomy = (SR*fMy_L - SL*fMy_R + SL*SR*(rhoR*vyR - rhoL*vyL)) * inv_SRSL
                    fener = (SR*fE_L  - SL*fE_R  + SL*SR*(enR      - enL     )) * inv_SRSL
                    fbyf  = (SR*fB_L  - SL*fB_R  + SL*SR*(ByR      - ByL     )) * inv_SRSL
                end if

                flux_Mass(i,j)   = fmass
                flux_Momx(i,j)   = fmomx
                flux_Momy(i,j)   = fmomy
                flux_Energy(i,j) = fener
                flux_By(i,j)     = fbyf

            end do
        end do

    end subroutine hlld_flux

end module mhd_flux
