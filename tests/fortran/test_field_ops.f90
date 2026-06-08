program test_field_ops
    !
    ! Unit tests for mhd_field_ops (compute_curl_2d, compute_divB,
    ! compute_gradients from mhd_derivatives).
    !
    ! Tests:
    !   1. div(curl(Az)) = 0 exactly by discrete algebra (periodic BCs)
    !   2. Gradient of a linear field is exact to machine precision
    !   3. Gradient of a sinusoidal field matches analytical result to O(dx^2)
    !
    use mhd_config
    use mhd_field_ops
    use mhd_derivatives
    implicit none

    integer, parameter :: nx = 32
    integer, parameter :: ny = 32
    real(8), parameter :: tol_exact  = 1.0d-13   ! for algebraically exact results
    real(8), parameter :: tol_second = 0.1d0      ! for O(dx^2) truncation error
    ! tol_second rationale: for f=sin(2*pi*x) the central-difference gradient
    ! error is (2*pi)^3 * dx^2 / 6 ~ 0.04 at N=32. 0.1 gives comfortable margin.

    real(8) :: dx_t, x(nx), y(ny)
    real(8) :: Az(nx,ny), b_x(nx,ny), b_y(nx,ny), divB(nx,ny)
    real(8) :: f(nx,ny), f_pad(nx+2,ny+2), f_dx(nx,ny), f_dy(nx,ny)
    real(8) :: f_dx_exact(nx,ny), f_dy_exact(nx,ny)
    integer :: i, j, n_pass, n_fail

    n_pass = 0
    n_fail = 0
    N = nx

    ! Set periodic BCs for all tests
    BC_xlo = BC_PERIODIC;  BC_xhi = BC_PERIODIC
    BC_ylo = BC_PERIODIC;  BC_yhi = BC_PERIODIC

    dx_t = boxsize / dble(nx)
    do i = 1, nx
        x(i) = (i - 0.5d0) * dx_t
    end do
    do j = 1, ny
        y(j) = (j - 0.5d0) * dx_t
    end do

    ! -------------------------------------------------------------------
    ! Test 1: div(curl(Az)) = 0 exactly
    ! The discrete backward-difference curl followed by the discrete
    ! backward-difference divergence cancels exactly by algebra,
    ! regardless of the field values. This confirms the CT scheme is
    ! algebraically consistent.
    ! -------------------------------------------------------------------
    do j = 1, ny
        do i = 1, nx
            Az(i,j) = sin(twoPi*x(i)) * cos(twoPi*y(j)) &
                    + 0.3d0*cos(fourPi*x(i)) * sin(twoPi*y(j))
        end do
    end do

    call compute_curl_2d(Az, dx_t, nx, ny, b_x, b_y)
    call compute_divB(b_x, b_y, dx_t, nx, ny, divB)

    call check("div(curl(Az)) = 0 (interior)", &
               maxval(abs(divB(2:nx, 2:ny))) < tol_exact, n_pass, n_fail)
    call check("div(curl(Az)) = 0 (boundaries)", &
               max(maxval(abs(divB(1,:))), maxval(abs(divB(:,1)))) < tol_exact, &
               n_pass, n_fail)

    ! -------------------------------------------------------------------
    ! Test 2: Gradient of a linear field is exact
    ! f(x,y) = 3*x + 2*y  =>  df/dx = 3, df/dy = 2 everywhere
    ! Central differences are exact for linear functions regardless of dx.
    ! -------------------------------------------------------------------
    do j = 1, ny
        do i = 1, nx
            f(i,j) = 3.0d0*x(i) + 2.0d0*y(j)
        end do
    end do

    ! Fill padded array with periodic ghost cells manually
    f_pad(2:nx+1, 2:ny+1) = f
    f_pad(1,    2:ny+1) = f(nx, :)   ! left ghost
    f_pad(nx+2, 2:ny+1) = f(1,  :)   ! right ghost
    f_pad(2:nx+1, 1)    = f(:, ny)   ! bottom ghost
    f_pad(2:nx+1, ny+2) = f(:, 1)    ! top ghost
    f_pad(1,    1)    = f(nx, ny);  f_pad(nx+2, 1)    = f(1, ny)
    f_pad(1,    ny+2) = f(nx, 1);   f_pad(nx+2, ny+2) = f(1, 1)

    call compute_gradients(f_pad, dx_t, nx, ny, f_dx, f_dy)

    ! Check only interior cells (i=2..nx-1, j=2..ny-1).
    ! At boundary cells (i=1 and i=nx), the periodic ghost wraps f(nx,:)
    ! to the left, which is wrong for a non-periodic linear field -- this
    ! is correct ghost-cell behaviour, not a bug. Central differences at
    ! interior cells are algebraically exact for linear functions.
    call check("linear grad: df/dx = 3", &
               maxval(abs(f_dx(2:nx-1, 2:ny-1) - 3.0d0)) < tol_exact, n_pass, n_fail)
    call check("linear grad: df/dy = 2", &
               maxval(abs(f_dy(2:nx-1, 2:ny-1) - 2.0d0)) < tol_exact, n_pass, n_fail)

    ! -------------------------------------------------------------------
    ! Test 3: Gradient of a sinusoidal field matches analytical O(dx^2)
    ! f(x,y) = sin(2*pi*x)*cos(2*pi*y)
    ! df/dx = 2*pi*cos(2*pi*x)*cos(2*pi*y)
    ! df/dy = -2*pi*sin(2*pi*x)*sin(2*pi*y)
    ! -------------------------------------------------------------------
    do j = 1, ny
        do i = 1, nx
            f(i,j)        = sin(twoPi*x(i)) * cos(twoPi*y(j))
            f_dx_exact(i,j) = twoPi * cos(twoPi*x(i)) * cos(twoPi*y(j))
            f_dy_exact(i,j) = -twoPi * sin(twoPi*x(i)) * sin(twoPi*y(j))
        end do
    end do

    f_pad(2:nx+1, 2:ny+1) = f
    f_pad(1,    2:ny+1) = f(nx, :)
    f_pad(nx+2, 2:ny+1) = f(1,  :)
    f_pad(2:nx+1, 1)    = f(:, ny)
    f_pad(2:nx+1, ny+2) = f(:, 1)
    f_pad(1,    1)    = f(nx, ny);  f_pad(nx+2, 1)    = f(1, ny)
    f_pad(1,    ny+2) = f(nx, 1);   f_pad(nx+2, ny+2) = f(1, 1)

    call compute_gradients(f_pad, dx_t, nx, ny, f_dx, f_dy)

    call check("sinusoidal grad: df/dx (O(dx^2))", &
               maxval(abs(f_dx - f_dx_exact)) < tol_second, n_pass, n_fail)
    call check("sinusoidal grad: df/dy (O(dx^2))", &
               maxval(abs(f_dy - f_dy_exact)) < tol_second, n_pass, n_fail)

    ! -------------------------------------------------------------------
    ! Summary
    ! -------------------------------------------------------------------
    write(*, '(A, I0, A, I0, A)') &
        "field_ops: ", n_pass, " passed, ", n_fail, " failed"
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

end program test_field_ops
