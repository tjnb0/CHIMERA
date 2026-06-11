module mhd_change_states
    use mhd_config
    !----------------------------------------------------------------------------
    ! Purpose: Subroutines to converting between primitive and conserved variables
    !          (density, velocity, pressure,  magnetic field) <->
    !          (mass, momentum, total energy, magnetic field)
    !-----------------------------------------------------------------------------

    implicit none
    private
    public get_conserved, get_primitive, apply_conserved_floors
    
contains


    subroutine get_conserved(rho, vx, vy, P, Bx, By, gamma, vol, nx, ny, &
                             Mass, Momx, Momy, Energy)
    ! 
    !   Convert primitive variables to conserved variables.
    !
    !   Inputs:
    !       - rho : (nx, ny) array
    !       - vx  : (nx, ny) array
    !       - vy  : (nx, ny) array
    !       - P   : (nx, ny) array
    !       - Bx  : (nx, ny) array
    !       - By  : (nx, ny) array
    !       - gamma : ideal gas gamma (scalar)
    !       - vol   : cell volume (scalar)
    !
    !   Outputs:
    !       - Mass   : (nx, ny) array;  rho * vol
    !       - Momx   : (nx, ny) array;  rho * vx * vol
    !       - Momy   : (nx, ny) array;  rho * vy * vol
    !       - Energy : (nx, ny) array;  (internal + kinetic + magnetic) * vol
    !
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: rho(nx, ny), vx(nx, ny), vy(nx, ny)
        real(8), intent(in) :: P(nx, ny), Bx(nx, ny), By(nx, ny)
        real(8), intent(in) :: gamma, vol
        real(8), intent(out) :: Mass(nx,ny), Momx(nx,ny), Momy(nx,ny), Energy(nx,ny)
        real(8) :: halfB2(nx, ny), rho_cp(nx, ny), P_cp(nx, ny)

        ! For realism
        rho_cp = max(rho, rho_floor)
        P_cp   = max(P,   P_floor)

        Mass   = rho_cp  * vol
        Momx   = Mass * vx
        Momy   = Mass * vy
        halfB2 = 0.5d0 * (Bx*Bx + By*By)          ! Magnetic: 0.5 * (Bx^2 + By^2)

        ! Total energy
        Energy = ((P_cp - halfB2) / (gamma - 1.d0) & ! Internal: (P - 0.5*B^2) / (gamma-1)
                + 0.5d0 * rho_cp * (vx*vx + vy*vy) & ! Kinetic : 0.5 * rho * (vx^2 + vy^2)
                + halfB2) * vol       

    end subroutine get_conserved


    subroutine get_primitive(Mass, Momx, Momy, Energy, Bx, By, gamma, vol, nx, ny, &
                             rho, vx, vy, P)
    ! 
    !   Convert conserved variables to primitive variables.
    !
    !   Inputs:
    !       - Mass   : (nx, ny) array
    !       - Momx   : (nx, ny) array
    !       - Momy   : (nx, ny) array
    !       - Energy : (nx, ny) array
    !       - Bx     : (nx, ny) array
    !       - By     : (nx, ny) array
    !       - gamma  : scalar
    !       - vol    : scalar
    !
    !   Outputs:
    !       - rho : (nx, ny) array;  Mass / vol
    !       - vx  : (nx, ny) array;  Momx / (rho * vol)
    !       - vy  : (nx, ny) array;  Momy / (rho * vol)
    !       - P   : (nx, ny) array;  Pressure
    ! 

        integer, intent(in) :: nx, ny
        real(8), intent(in) :: Mass(nx, ny), Momx(nx, ny), Momy(nx, ny)
        real(8), intent(in) :: Energy(nx, ny), Bx(nx, ny), By(nx, ny)
        real(8), intent(in) :: gamma, vol
        real(8), intent(out) :: rho(nx, ny), vx(nx, ny), vy(nx, ny), P(nx, ny)
        real(8) :: halfB2(nx, ny)


        rho = max(Mass / vol, rho_floor) ! apply safety check
        vx  = Momx / (rho * vol)         ! use safe rho (not Mass) to get vel
        vy  = Momy / (rho * vol)         ! use safe rho (not Mass) to get vel
        halfB2 = 0.5d0 * (Bx*Bx + By*By)

        !   P = (Energy/vol - K.E. - Mag. Energy) * (gamma-1) + Mag. Pressure
        !     - Kinetic       = 0.5 * rho * (vx^2 + vy^2)
        !     - Mag. Energy   = 0.5 * (Bx^2 + By^2)
        !     - Mag. Pressure = 0.5 * (Bx^2 + By^2)
        P = (Energy / vol - 0.5d0*rho*(vx*vx + vy*vy) - halfB2) &
          * (gamma - 1.d0) + halfB2

        ! Issue 2 fix: absolute floor first, then proportional floor.
        ! At high Mach / low beta, thermal pressure is a small residual of
        ! large numbers. The proportional floor ensures p is always a
        ! meaningful fraction of the local kinetic + magnetic energy,
        ! preventing cancellation errors from corrupting the solution.
        ! e_floor_frac is tunable in mhd_config.f90.
        P = max(P, P_floor)
        P = max(P, halfB2 + e_floor_frac * (0.5d0*rho*(vx*vx + vy*vy) + halfB2))
        

    end subroutine get_primitive


    subroutine apply_conserved_floors(Mass, Momx, Momy, Energy, Bx, By, vol, nx, ny)
    !
    !   Enforce physical limits on conserved variables after update_conserved.
    !   Called once per timestep before get_primitive.
    !
    !   Checks applied:
    !       1. Mass floor: negative mass is unphysical; zero momentum and
    !          reset energy to magnetic-only floor when mass is floored.
    !       2. Energy floor: energy below the magnetic floor implies negative
    !          thermal pressure regardless of momentum; clamp it.
    !       3. Momentum limiting: cap implied velocity at V_MAX to prevent
    !          near-vacuum cells with residual momentum from blowing up the CFL
    !          and the reconstruction prediction step.
    !
    !
        integer, intent(in)    :: nx, ny
        real(8), intent(inout) :: Mass(nx,ny), Momx(nx,ny), Momy(nx,ny)
        real(8), intent(inout) :: Energy(nx,ny)
        real(8), intent(in)    :: Bx(nx,ny), By(nx,ny)
        real(8), intent(in)    :: vol

        real(8) :: halfB2(nx,ny), E_mag_floor(nx,ny)
        real(8) :: v_mag(nx,ny), scale(nx,ny)

        real(8), parameter :: V_MAX = 50.0d0   ! velocity cap for near-vacuum cells

        ! Magnetic-only energy floor: E >= 0.5*(Bx^2+By^2)*vol
        halfB2      = 0.5d0*(Bx*Bx + By*By)
        E_mag_floor = halfB2 * vol

        ! Check mass floor 
        ! Negative mass -> velocity and thermal pressure blow up. Zero momentum
        ! and sync energy to a magnetic state
        where (Mass < rho_floor * vol)
            Momx   = 0.0d0
            Momy   = 0.0d0
            Mass   = rho_floor * vol
            Energy = E_mag_floor
        end where

        ! Check energy floor
        ! Energy below the magnetic floor means negative thermal pressure
        where (Energy < E_mag_floor)
            Energy = E_mag_floor
        end where

        ! Check momentum limiting
        ! If v = |Mom| / Mass > V_MAX, the cell is near-vacuum and gives
        ! non-physical flow. Scale components and preserve direction.
        v_mag = sqrt(Momx*Momx + Momy*Momy) / Mass
        where (v_mag > V_MAX)
            scale = V_MAX * Mass / sqrt(Momx*Momx + Momy*Momy)
            Momx  = Momx * scale
            Momy  = Momy * scale
        end where

    end subroutine apply_conserved_floors
    
end module mhd_change_states