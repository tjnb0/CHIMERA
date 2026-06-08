program test_change_states
    !
    ! Unit tests for mhd_change_states (get_conserved / get_primitive).
    !
    ! Tests:
    !   1. Round-trip: primitives -> conserved -> primitives recovers original values
    !   2. Density floor: low-density cells do not produce velocity blow-up
    !   3. Proportional pressure floor: activates when thermal pressure << kinetic energy
    !
    use mhd_config
    use mhd_change_states
    implicit none

    integer,  parameter :: nx    = 4
    integer,  parameter :: ny    = 4
    real(8),  parameter :: tol   = 1.0d-10
    real(8),  parameter :: gamma_t = 5.0d0 / 3.0d0
    real(8),  parameter :: vol_t   = (1.0d0/64.0d0)**2  ! dx^2 for N=64

    real(8) :: rho_in(nx,ny), vx_in(nx,ny), vy_in(nx,ny)
    real(8) :: P_in(nx,ny),   Bx_in(nx,ny), By_in(nx,ny)
    real(8) :: Mass(nx,ny), Momx(nx,ny), Momy(nx,ny), Energy(nx,ny)
    real(8) :: rho_out(nx,ny), vx_out(nx,ny), vy_out(nx,ny), P_out(nx,ny)
    real(8) :: halfB2_val, KE_val, P_floor_expected

    integer :: n_pass, n_fail

    n_pass = 0
    n_fail = 0
    N = nx  ! set the module variable

    ! -------------------------------------------------------------------
    ! Test 1: Round-trip for a typical plasma state
    ! get_conserved followed by get_primitive should recover the input
    ! -------------------------------------------------------------------
    rho_in = 1.5d0
    vx_in  = 0.3d0
    vy_in  = -0.2d0
    Bx_in  = 0.5d0
    By_in  = 0.3d0
    ! Set total pressure: P* = p_thermal + halfB2
    P_in   = 1.0d0 + 0.5d0*(Bx_in(1,1)**2 + By_in(1,1)**2)

    call get_conserved(rho_in, vx_in, vy_in, P_in, Bx_in, By_in, &
                       gamma_t, vol_t, nx, ny, Mass, Momx, Momy, Energy)
    call get_primitive(Mass, Momx, Momy, Energy, Bx_in, By_in, &
                       gamma_t, vol_t, nx, ny, rho_out, vx_out, vy_out, P_out)

    call check("round-trip rho", maxval(abs(rho_out - rho_in)) < tol,  n_pass, n_fail)
    call check("round-trip vx",  maxval(abs(vx_out  - vx_in))  < tol,  n_pass, n_fail)
    call check("round-trip vy",  maxval(abs(vy_out  - vy_in))  < tol,  n_pass, n_fail)
    call check("round-trip P",   maxval(abs(P_out   - P_in))   < tol,  n_pass, n_fail)

    ! -------------------------------------------------------------------
    ! Test 2: Density floor prevents velocity blow-up (Issue 1 fix)
    ! Before the fix: vx = Momx/Mass used the un-floored Mass, giving
    ! enormous velocities when rho < rho_floor.
    ! After the fix:  vx = Momx/(rho*vol) uses the floored rho.
    ! -------------------------------------------------------------------
    rho_in = 1.0d-20   ! far below rho_floor = 1e-12
    vx_in  = 1.0d0
    vy_in  = 0.0d0
    P_in   = 1.0d0
    Bx_in  = 0.0d0
    By_in  = 0.0d0

    call get_conserved(rho_in, vx_in, vy_in, P_in, Bx_in, By_in, &
                       gamma_t, vol_t, nx, ny, Mass, Momx, Momy, Energy)
    call get_primitive(Mass, Momx, Momy, Energy, Bx_in, By_in, &
                       gamma_t, vol_t, nx, ny, rho_out, vx_out, vy_out, P_out)

    call check("density floor applied",       minval(rho_out) >= rho_floor, n_pass, n_fail)
    call check("velocity finite after floor", maxval(abs(vx_out)) < 1.0d6,  n_pass, n_fail)

    ! -------------------------------------------------------------------
    ! Test 3: Proportional pressure floor activates at high kinetic energy
    ! (Issue 2 fix).
    ! Construct a state where thermal pressure is zero (all energy is
    ! kinetic + magnetic). After get_primitive, P should be above the
    ! proportional floor: P >= halfB2 + e_floor_frac*(KE + halfB2).
    ! -------------------------------------------------------------------
    rho_in = 1.0d0
    vx_in  = 1.0d3   ! very fast: KE >> thermal energy
    vy_in  = 0.0d0
    Bx_in  = 0.5d0
    By_in  = 0.0d0

    halfB2_val = 0.5d0 * (Bx_in(1,1)**2 + By_in(1,1)**2)
    KE_val     = 0.5d0 * rho_in(1,1) * (vx_in(1,1)**2 + vy_in(1,1)**2)

    ! Build conserved state with zero thermal energy
    Mass   = rho_in * vol_t
    Momx   = rho_in * vx_in * vol_t
    Momy   = 0.0d0
    Energy = (KE_val + halfB2_val) * vol_t   ! no thermal contribution

    call get_primitive(Mass, Momx, Momy, Energy, Bx_in, By_in, &
                       gamma_t, vol_t, nx, ny, rho_out, vx_out, vy_out, P_out)

    P_floor_expected = halfB2_val + e_floor_frac * (KE_val + halfB2_val)

    call check("proportional floor: P > 0",          minval(P_out) > 0.0d0,              n_pass, n_fail)
    call check("proportional floor: P >= floor",     minval(P_out) >= P_floor_expected * (1.0d0 - 1.0d-10), n_pass, n_fail)

    ! -------------------------------------------------------------------
    ! Summary
    ! -------------------------------------------------------------------
    write(*, '(A, I0, A, I0, A)') &
        "change_states: ", n_pass, " passed, ", n_fail, " failed"
    if (n_fail > 0) stop 1

contains

    subroutine check(name, passed, n_pass, n_fail)
        character(len=*), intent(in)    :: name
        logical,          intent(in)    :: passed
        integer,          intent(inout) :: n_pass, n_fail
        if (passed) then
            write(*, '(A, A)') "  PASS  ", name
            n_pass = n_pass + 1
        else
            write(*, '(A, A)') "  FAIL  ", name
            n_fail = n_fail + 1
        end if
    end subroutine check

end program test_change_states
