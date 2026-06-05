module mhd_derivatives
    use mhd_config
    !---------------------------------------------------------------------------
    ! Purpose: Subroutines for computing spatial derivatives, applying slope
    !          limiters, and performing spatial extrapolation to cell faces for
    !          finite-volume MHD.
    !---------------------------------------------------------------------------

    implicit none
    private
    public compute_gradients, apply_slope_limiter, reconstruction
    
contains

    subroutine compute_gradients(f, dx, nx, ny, f_dx, f_dy)
    !
    !   Calculate the gradients of a field with periodic boundaries.
    !
    !   Inputs:
    !       - f  : (nx, ny) array;  field
    !       - dx : scalar;          cell size
    !
    !   Outputs:
    !       - f_dx : (nx, ny) array;  df/dx
    !       - f_dy : (nx, ny) array;  df/dy
    !
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: f(nx, ny)
        real(8), intent(in) :: dx
        real(8), intent(out) :: f_dx(nx, ny), f_dy(nx, ny)
        real(8) :: inv_2dx

        inv_2dx = 1.0d0 / (2.0d0 * dx)

        ! x-derivative (df/dx) 
        f_dx(2:nx-1,:) = (f(3:nx,:) - f(1:nx-2,:)) * inv_2dx ! Interior points (cent. diff)
        f_dx(1,:)      = (f(2,:)    - f(nx,:))     * inv_2dx ! Left boundary   (periodic)
        f_dx(nx,:)     = (f(1,:)    - f(nx-1,:))   * inv_2dx ! Right boundary  (periodic)

        ! y-derivative (df/dy) 
        f_dy(:,2:ny-1) = (f(:,3:ny) - f(:,1:ny-2)) * inv_2dx ! Interior points (cent. diff)
        f_dy(:,1)      = (f(:,2)    - f(:,ny))     * inv_2dx ! Lower boundary  (periodic)
        f_dy(:,ny)     = (f(:,1)    - f(:,ny-1))   * inv_2dx ! Upper boundary  (periodic)

    end subroutine compute_gradients


    subroutine apply_slope_limiter(f, dx, nx, ny, f_dx, f_dy)
    ! 
    !   Apply slope limiter to x- and y- directional slopes
    !
    !   Inputs:
    !       - f  : scalar field (nx, ny)
    !       - dx : grid spacing
    !       - nx : x grid size
    !       - ny : y grid size
    !
    ! InOut:
    !       - f_dx : slope to be limited (nx)
    !       - f_dy : slope to be limited (ny)
    !
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: dx
        real(8), intent(in) :: f(nx,ny)
        real(8), intent(out) :: f_dx(nx,ny), f_dy(nx,ny)
        integer :: i, j, ip1, im1, jp1, jm1
        real(8) :: dfR, dfL, dfRdfL, inv_dx, den

        inv_dx = 1.0d0 / dx

        !$omp parallel do private(i,j,ip1,im1,dfR,dfL,dfRdfL)
        do j = 1, ny
            do i = 1, nx
                ip1 = mod(i, nx) + 1
                im1 = mod(i-2+nx, nx) + 1

                dfL    = (f(i,j) - f(im1,j)) * inv_dx
                dfR    = (f(ip1,j) - f(i,j)) * inv_dx
                dfRdfL = dfL * dfR
                den    = dfL + dfR
                if (dfRdfL <= 0.d0 .or. abs(den) < 1d-12) then
                    f_dx(i,j) = 0.d0
                else
                    f_dx(i,j) = (2.d0 * dfRdfL)/(dfL + dfR)
                end if
            end do
        end do

        !$omp parallel do private(i,j,jp1,jm1,dfR,dfL,dfRdfL)
        do j = 1, ny
            do i = 1, nx
                jp1 = mod(j, ny) + 1
                jm1 = mod(j-2+ny, ny) + 1

                dfL = (f(i,j) - f(i,jm1)) * inv_dx
                dfR = (f(i,jp1) - f(i,j)) * inv_dx
                dfRdfL = dfL * dfR
                den    = dfL + dfR
                if (dfRdfL <= 0.d0 .or. abs(den) < 1d-12) then
                    f_dy(i,j) = 0.d0
                else
                    f_dy(i,j) = (2.d0 * dfRdfL)/(dfL + dfR)
                end if
            end do
        end do
    end subroutine apply_slope_limiter


    subroutine reconstruction(f, f_dx, f_dy, dx, nx, ny, f_XL, f_XR, f_YL, f_YR, &
                                   rho_or_p)
    ! 
    !   Performs MUSCL- or MOOD-scheme spatial extrapolation to cell faces using 
    !   slope-limited derivatives
    !
    !   Inputs:
    !       - f    : field (nx, ny)
    !       - f_dx : limited x-derivative of f, (nx, ny)
    !       - f_dy : limited y-derivative of f, (nx, ny)
    !       - dx   : cell size
    !       - nx   : x grid size
    !       - ny   : y grid size
    !
    !   Outputs:
    !       - f_XL : extrapolated values on left  face in x-direction
    !       - f_XR : extrapolated values on right face in x-direction
    !       - f_YL : extrapolated values on left  face in y-direction
    !       - f_YR : extrapolated values on right face in y-direction
    !
    !  Source:
    !       - https://en.wikipedia.org/wiki/MUSCL_scheme
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

        
        ! 1. Use MUSCL reconstruction to make candidates
        f_XR   = f + 0.5d0 * dx * f_dx ! Right face
        f_YR   = f + 0.5d0 * dx * f_dy ! Top face
        temp_x = f - 0.5d0 * dx * f_dx ! Left face (before shift)
        temp_y = f - 0.5d0 * dx * f_dy ! Bottom face (before shift)
        f_XL = cshift(temp_x, +1, 1)   ! Periodic shift right
        f_YL = cshift(temp_y, +1, 2)   ! Periodic shift up

        ! 2. Optionally use MOOD reconstruction for better accuracy
        if (upgrade_2_MOOD) then
            !$OMP parallel do private(i,j,ip1,im1,jp1,jm1,localmin,localmax) &
            !$OMP shared(f, f_XL, f_XR, f_YL, f_YR, nx, ny, rho_or_p)        
            do j = 1, ny
                do i = 1, nx
                    ! Compute indices for periodic neighbors (von Neumann stencil)
                    ip1 = mod(i, nx) + 1      ! i+1
                    im1 = mod(i-2, nx) + 1    ! i-1
                    jp1 = mod(j, ny) + 1      ! j+1
                    jm1 = mod(j-2, ny) + 1    ! j-1

                    ! Local min/max over cell and neighbors (wider stencil)
                    localmin = min(f(i,j), f(ip1,j), f(im1,j), f(i,jp1), f(i,jm1))
                    localmax = max(f(i,j), f(ip1,j), f(im1,j), f(i,jp1), f(i,jm1))

                    ! Non-negative density/pressure
                    if (rho_or_p) then
                        localmin = max(0.1d-8, localmin)
                    end if

                    ! Check x-direction reconstructions
                    if (f_XL(i,j) < localmin .or. f_XL(i,j) > localmax) then
                        ! Fallback to first-order (Godunov-like)
                        f_XL(i,j) = f(i,j)
                    end if
                    if (f_XR(i,j) < localmin .or. f_XR(i,j) > localmax) then
                        f_XR(i,j) = f(i,j)
                    end if

                    ! Check y-direction reconstructions
                    if (f_YL(i,j) < localmin .or. f_YL(i,j) > localmax) then
                        f_YL(i,j) = f(i,j)
                    end if
                    if (f_YR(i,j) < localmin .or. f_YR(i,j) > localmax) then
                        f_YR(i,j) = f(i,j)
                    end if

                end do
            end do
            !$omp end parallel do
        end if

        ! For safer reconstruction
        eps = 1e-8
        if (rho_or_p) then
            f_XL = max(f_XL, eps)
            f_XR = max(f_XR, eps)
            f_YL = max(f_YL, eps)
            f_YR = max(f_YR, eps)
        end if

    end subroutine reconstruction


end module mhd_derivatives