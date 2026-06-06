module mhd_derivatives
    use mhd_config
    !---------------------------------------------------------------------------
    ! Purpose: Subroutines for computing spatial derivatives, applying slope
    !          limiters, and performing spatial extrapolation to cell faces for
    !          finite-volume MHD.
    !---------------------------------------------------------------------------

    implicit none
    private
    public compute_gradients, apply_slope_limiter, reconstruction, thermal_pressure_check
    
contains

    subroutine compute_gradients(f_pad, dx, nx, ny, f_dx, f_dy)
    !
    !   Compute cell-centred gradients using a ghost-cell padded field.
    !   Physical cell (i,j) maps to f_pad(i+1, j+1). Ghost cells in rows/
    !   columns 1 and nx+2/ny+2 encode the BC (periodic wrap, zero-gradient
    !   outflow, or prescribed inflow) and are filled by fill_ghost_cells()
    !   before this call. No special boundary cases are needed here.
    !
    !   Inputs:
    !       - f_pad : (nx+2, ny+2) padded field
    !       - dx    : cell size
    !
    !   Outputs:
    !       - f_dx  : (nx, ny) df/dx
    !       - f_dy  : (nx, ny) df/dy
    !
        integer, intent(in)  :: nx, ny
        real(8), intent(in)  :: f_pad(nx+2, ny+2)
        real(8), intent(in)  :: dx
        real(8), intent(out) :: f_dx(nx, ny), f_dy(nx, ny)
        real(8) :: inv_2dx

        inv_2dx = 1.0d0 / (2.0d0 * dx)

        ! Central differences over the full domain - ghost cells handle all BCs.
        ! Left neighbor of cell i   is f_pad(i,   j+1).
        ! Right neighbor of cell i  is f_pad(i+2, j+1).
        ! Bottom neighbor of cell j is f_pad(i+1, j  ).
        ! Top neighbor of cell j    is f_pad(i+1, j+2).
        f_dx = (f_pad(3:nx+2, 2:ny+1) - f_pad(1:nx,   2:ny+1)) * inv_2dx
        f_dy = (f_pad(2:nx+1, 3:ny+2) - f_pad(2:nx+1, 1:ny  )) * inv_2dx

    end subroutine compute_gradients


    subroutine apply_slope_limiter(f_pad, dx, nx, ny, f_dx, f_dy)
    !
    !   Apply Van Leer harmonic mean slope limiter using ghost-cell padded field.
    !   Ghost cells encode the BC so no special boundary handling is needed.
    !
    !   Inputs:
    !       - f_pad : (nx+2, ny+2) padded field
    !       - dx    : grid spacing
    !
    !   InOut:
    !       - f_dx  : (nx, ny) x-slope to be limited
    !       - f_dy  : (nx, ny) y-slope to be limited
    !
        integer, intent(in)  :: nx, ny
        real(8), intent(in)  :: dx
        real(8), intent(in)  :: f_pad(nx+2, ny+2)
        real(8), intent(out) :: f_dx(nx, ny), f_dy(nx, ny)
        integer :: i, j
        real(8) :: dfR, dfL, dfRdfL, inv_dx, den

        inv_dx = 1.0d0 / dx

        ! x-direction slopes.
        ! Physical cell (i,j) = f_pad(i+1, j+1).
        ! Left neighbor  = f_pad(i,   j+1).
        ! Right neighbor = f_pad(i+2, j+1).
        !$omp parallel do private(i,j,dfR,dfL,dfRdfL,den)
        do j = 1, ny
            do i = 1, nx
                dfL    = (f_pad(i+1,j+1) - f_pad(i,  j+1)) * inv_dx
                dfR    = (f_pad(i+2,j+1) - f_pad(i+1,j+1)) * inv_dx
                dfRdfL = dfL * dfR
                den    = dfL + dfR
                if (dfRdfL <= 0.d0 .or. abs(den) < 1d-12) then
                    f_dx(i,j) = 0.d0
                else
                    f_dx(i,j) = (2.d0 * dfRdfL) / den
                end if
            end do
        end do
        !$omp end parallel do

        ! y-direction slopes.
        ! Bottom neighbor = f_pad(i+1, j  ).
        ! Top neighbor    = f_pad(i+1, j+2).
        !$omp parallel do private(i,j,dfR,dfL,dfRdfL,den)
        do j = 1, ny
            do i = 1, nx
                dfL    = (f_pad(i+1,j+1) - f_pad(i+1,j  )) * inv_dx
                dfR    = (f_pad(i+1,j+2) - f_pad(i+1,j+1)) * inv_dx
                dfRdfL = dfL * dfR
                den    = dfL + dfR
                if (dfRdfL <= 0.d0 .or. abs(den) < 1d-12) then
                    f_dy(i,j) = 0.d0
                else
                    f_dy(i,j) = (2.d0 * dfRdfL) / den
                end if
            end do
        end do
        !$omp end parallel do

    end subroutine apply_slope_limiter


    subroutine reconstruction(f, f_dx, f_dy, dx, nx, ny, f_XL, f_XR, f_YL, f_YR, &
                                   rho_or_p)
    !
    !   MUSCL/MOOD spatial extrapolation to cell faces.
    !   cshift calls replaced with explicit BC-aware boundary handling
    !   using the per-side flags from mhd_config.
    !
    !   Inputs:
    !       - f    : field (nx, ny)
    !       - f_dx : limited x-derivative of f, (nx, ny)
    !       - f_dy : limited y-derivative of f, (nx, ny)
    !       - dx   : cell size
    !
    !   Outputs:
    !       - f_XL : right state at the right x-face of each cell
    !       - f_XR : left  state at the right x-face of each cell
    !       - f_YL : right state at the top   y-face of each cell
    !       - f_YR : left  state at the top   y-face of each cell
    !
        integer, intent(in)  :: nx, ny
        real(8), intent(in)  :: dx
        real(8), intent(in)  :: f(nx, ny), f_dx(nx, ny), f_dy(nx, ny)
        real(8), intent(out) :: f_XL(nx, ny), f_XR(nx, ny)
        real(8), intent(out) :: f_YL(nx, ny), f_YR(nx, ny)
        logical, intent(in)  :: rho_or_p

        real(8) :: temp_x(nx, ny), temp_y(nx, ny)
        real(8) :: localmin, localmax, eps
        integer :: i, j, ip1, im1, jp1, jm1

        ! 1. MUSCL reconstruction candidates
        f_XR   = f + 0.5d0 * dx * f_dx   ! left  state at right x-face
        f_YR   = f + 0.5d0 * dx * f_dy   ! left  state at top   y-face
        temp_x = f - 0.5d0 * dx * f_dx   ! right state at right x-face (pre-shift)
        temp_y = f - 0.5d0 * dx * f_dy   ! right state at top   y-face (pre-shift)

        ! X right state: interior uses right neighbor; boundary uses BC.
        f_XL(1:nx-1, :) = temp_x(2:nx, :)
        select case (bc_xhi)
            case (BC_PERIODIC)
                f_XL(nx, :) = temp_x(1, :)   ! wrap
            case default                       ! outflow/fixed/inflow: zero-gradient
                f_XL(nx, :) = f(nx, :)
        end select

        ! Y right state: interior uses upper neighbor; boundary uses BC.
        f_YL(:, 1:ny-1) = temp_y(:, 2:ny)
        select case (bc_yhi)
            case (BC_PERIODIC)
                f_YL(:, ny) = temp_y(:, 1)   ! wrap
            case default
                f_YL(:, ny) = f(:, ny)
        end select

        ! 2. Optional MOOD fallback to first order at troubled cells
        if (upgrade_2_MOOD) then
            !$OMP parallel do private(i,j,ip1,im1,jp1,jm1,localmin,localmax) &
            !$OMP shared(f, f_XL, f_XR, f_YL, f_YR, nx, ny, rho_or_p)
            do j = 1, ny
                do i = 1, nx
                    ! BC-aware neighbor indices for MOOD stencil
                    select case (bc_xhi)
                        case (BC_PERIODIC); ip1 = mod(i, nx) + 1
                        case default;       ip1 = min(i + 1, nx)
                    end select
                    select case (bc_xlo)
                        case (BC_PERIODIC); im1 = mod(i-2+nx, nx) + 1
                        case default;       im1 = max(i - 1, 1)
                    end select
                    select case (bc_yhi)
                        case (BC_PERIODIC); jp1 = mod(j, ny) + 1
                        case default;       jp1 = min(j + 1, ny)
                    end select
                    select case (bc_ylo)
                        case (BC_PERIODIC); jm1 = mod(j-2+ny, ny) + 1
                        case default;       jm1 = max(j - 1, 1)
                    end select

                    localmin = min(f(i,j), f(ip1,j), f(im1,j), f(i,jp1), f(i,jm1))
                    localmax = max(f(i,j), f(ip1,j), f(im1,j), f(i,jp1), f(i,jm1))

                    if (rho_or_p) localmin = max(0.1d-8, localmin)

                    if (f_XL(i,j) < localmin .or. f_XL(i,j) > localmax) f_XL(i,j) = f(i,j)
                    if (f_XR(i,j) < localmin .or. f_XR(i,j) > localmax) f_XR(i,j) = f(i,j)
                    if (f_YL(i,j) < localmin .or. f_YL(i,j) > localmax) f_YL(i,j) = f(i,j)
                    if (f_YR(i,j) < localmin .or. f_YR(i,j) > localmax) f_YR(i,j) = f(i,j)
                end do
            end do
            !$omp end parallel do
        end if

        ! Floor for density and pressure
        eps = 1.0d-8
        if (rho_or_p) then
            f_XL = max(f_XL, eps);  f_XR = max(f_XR, eps)
            f_YL = max(f_YL, eps);  f_YR = max(f_YR, eps)
        end if

    end subroutine reconstruction

    subroutine thermal_pressure_check(nx, ny,                        &
                                       P_XL,  P_XR,  P_YL,  P_YR,   &
                                       Bx_XL, Bx_XR, Bx_YL, Bx_YR, &
                                       By_XL, By_XR, By_YL, By_YR, &
                                       rho_XL, rho_XR, rho_YL, rho_YR, &
                                       vx_XL,  vx_XR,  vx_YL,  vx_YR,  &
                                       vy_XL,  vy_XR,  vy_YL,  vy_YR,  &
                                       rho_0, vx_0, vy_0, P_0, Bx_0, By_0)
    !
    !   Thermal-pressure MOOD criterion.
    !
    !   For any reconstructed face where the implied thermal pressure
    !   p = P_face - 0.5*(Bx_face^2 + By_face^2) drops below P_floor,
    !   all six fields at that face are reset to the predicted cell-centre
    !   values (first-order fallback). This extends the stencil-bounds
    !   check in reconstruction() to guard physical positivity.
    !
        integer, intent(in) :: nx, ny
        real(8), intent(inout) :: P_XL(nx,ny),   P_XR(nx,ny),   P_YL(nx,ny),   P_YR(nx,ny)
        real(8), intent(inout) :: Bx_XL(nx,ny),  Bx_XR(nx,ny),  Bx_YL(nx,ny),  Bx_YR(nx,ny)
        real(8), intent(inout) :: By_XL(nx,ny),  By_XR(nx,ny),  By_YL(nx,ny),  By_YR(nx,ny)
        real(8), intent(inout) :: rho_XL(nx,ny), rho_XR(nx,ny), rho_YL(nx,ny), rho_YR(nx,ny)
        real(8), intent(inout) :: vx_XL(nx,ny),  vx_XR(nx,ny),  vx_YL(nx,ny),  vx_YR(nx,ny)
        real(8), intent(inout) :: vy_XL(nx,ny),  vy_XR(nx,ny),  vy_YL(nx,ny),  vy_YR(nx,ny)
        real(8), intent(in)    :: rho_0(nx,ny), vx_0(nx,ny), vy_0(nx,ny)
        real(8), intent(in)    :: P_0(nx,ny),   Bx_0(nx,ny), By_0(nx,ny)

        real(8) :: p_th(nx,ny)

        ! X-left faces
        p_th = P_XL - 0.5d0*(Bx_XL**2 + By_XL**2)
        where (p_th < P_floor)
            P_XL  = P_0;  rho_XL = rho_0
            vx_XL = vx_0; vy_XL  = vy_0
            Bx_XL = Bx_0; By_XL  = By_0
        end where

        ! X-right faces
        p_th = P_XR - 0.5d0*(Bx_XR**2 + By_XR**2)
        where (p_th < P_floor)
            P_XR  = P_0;  rho_XR = rho_0
            vx_XR = vx_0; vy_XR  = vy_0
            Bx_XR = Bx_0; By_XR  = By_0
        end where

        ! Y-left faces
        p_th = P_YL - 0.5d0*(Bx_YL**2 + By_YL**2)
        where (p_th < P_floor)
            P_YL  = P_0;  rho_YL = rho_0
            vx_YL = vx_0; vy_YL  = vy_0
            Bx_YL = Bx_0; By_YL  = By_0
        end where

        ! Y-right faces
        p_th = P_YR - 0.5d0*(Bx_YR**2 + By_YR**2)
        where (p_th < P_floor)
            P_YR  = P_0;  rho_YR = rho_0
            vx_YR = vx_0; vy_YR  = vy_0
            Bx_YR = Bx_0; By_YR  = By_0
        end where

    end subroutine thermal_pressure_check


end module mhd_derivatives