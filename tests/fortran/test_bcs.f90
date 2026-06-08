program test_bcs
    !
    ! Unit tests for mhd_bc (fill_ghost_cells).
    !
    ! Tests:
    !   1. BC_PERIODIC: ghost cells match the opposite boundary
    !   2. BC_OUTFLOW:  ghost cells match the nearest interior cell
    !   3. Mixed BCs:   per-side flags applied independently
    !   4. Corners:     correct diagonal wrapping for periodic
    !
    use mhd_config
    use mhd_bc
    implicit none

    integer, parameter :: nx = 8
    integer, parameter :: ny = 8
    real(8), parameter :: tol = 1.0d-14

    real(8) :: f(nx, ny), f_pad(nx+2, ny+2)
    integer :: i, j, n_pass, n_fail

    n_pass = 0
    n_fail = 0
    N = nx

    ! Fill f with a recognisable non-uniform pattern
    do j = 1, ny
        do i = 1, nx
            f(i,j) = dble(i) + 10.0d0*dble(j)
        end do
    end do

    ! -------------------------------------------------------------------
    ! Test 1: BC_PERIODIC - ghost cells should match opposite boundary
    ! -------------------------------------------------------------------
    BC_xlo = BC_PERIODIC;  BC_xhi = BC_PERIODIC
    BC_ylo = BC_PERIODIC;  BC_yhi = BC_PERIODIC

    call fill_ghost_cells(f, f_pad)

    ! Left ghost should match right interior column
    call check("periodic: left ghost = right col",   &
               maxval(abs(f_pad(1, 2:ny+1) - f(nx, 1:ny))) < tol, n_pass, n_fail)
    ! Right ghost should match left interior column
    call check("periodic: right ghost = left col",   &
               maxval(abs(f_pad(nx+2, 2:ny+1) - f(1, 1:ny))) < tol, n_pass, n_fail)
    ! Bottom ghost should match top interior row
    call check("periodic: bottom ghost = top row",   &
               maxval(abs(f_pad(2:nx+1, 1) - f(1:nx, ny))) < tol, n_pass, n_fail)
    ! Top ghost should match bottom interior row
    call check("periodic: top ghost = bottom row",   &
               maxval(abs(f_pad(2:nx+1, ny+2) - f(1:nx, 1))) < tol, n_pass, n_fail)
    ! Interior should be unchanged
    call check("periodic: interior unchanged",       &
               maxval(abs(f_pad(2:nx+1, 2:ny+1) - f(1:nx, 1:ny))) < tol, n_pass, n_fail)
    ! Corner: lower-left should be upper-right physical cell (full diagonal wrap)
    call check("periodic: lower-left corner",        &
               abs(f_pad(1, 1) - f(nx, ny)) < tol, n_pass, n_fail)

    ! -------------------------------------------------------------------
    ! Test 2: BC_OUTFLOW - ghost cells should match nearest interior cell
    ! -------------------------------------------------------------------
    BC_xlo = BC_OUTFLOW;  BC_xhi = BC_OUTFLOW
    BC_ylo = BC_OUTFLOW;  BC_yhi = BC_OUTFLOW

    call fill_ghost_cells(f, f_pad)

    call check("outflow: left ghost = left col",     &
               maxval(abs(f_pad(1, 2:ny+1) - f(1, 1:ny))) < tol, n_pass, n_fail)
    call check("outflow: right ghost = right col",   &
               maxval(abs(f_pad(nx+2, 2:ny+1) - f(nx, 1:ny))) < tol, n_pass, n_fail)
    call check("outflow: bottom ghost = bottom row", &
               maxval(abs(f_pad(2:nx+1, 1) - f(1:nx, 1))) < tol, n_pass, n_fail)
    call check("outflow: top ghost = top row",       &
               maxval(abs(f_pad(2:nx+1, ny+2) - f(1:nx, ny))) < tol, n_pass, n_fail)
    ! Corner: nearest interior corner cell
    call check("outflow: lower-left corner = f(1,1)",   &
               abs(f_pad(1, 1) - f(1, 1)) < tol, n_pass, n_fail)
    call check("outflow: upper-right corner = f(N,N)",  &
               abs(f_pad(nx+2, ny+2) - f(nx, ny)) < tol, n_pass, n_fail)

    ! -------------------------------------------------------------------
    ! Test 3: Mixed BCs - x periodic, y outflow
    ! -------------------------------------------------------------------
    BC_xlo = BC_PERIODIC;  BC_xhi = BC_PERIODIC
    BC_ylo = BC_OUTFLOW;   BC_yhi = BC_OUTFLOW

    call fill_ghost_cells(f, f_pad)

    call check("mixed: left ghost periodic",         &
               maxval(abs(f_pad(1, 2:ny+1) - f(nx, 1:ny))) < tol, n_pass, n_fail)
    call check("mixed: bottom ghost outflow",        &
               maxval(abs(f_pad(2:nx+1, 1) - f(1:nx, 1))) < tol, n_pass, n_fail)
    ! Corner: x wraps, y doesn't
    call check("mixed: lower-left corner (x wraps)", &
               abs(f_pad(1, 1) - f(nx, 1)) < tol, n_pass, n_fail)

    ! -------------------------------------------------------------------
    ! Summary
    ! -------------------------------------------------------------------
    write(*, '(A, I0, A, I0, A)') &
        "bcs: ", n_pass, " passed, ", n_fail, " failed"
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

end program test_bcs
