module mhd_field_ops

    !---------------------------------------------------------------------------
    ! Purpose: Field operations for MHD solver: discrete curl, divergence, and 
    !          cell-averaging of vector fields. Used to update B-fields
    !---------------------------------------------------------------------------

    implicit none
    private
    public compute_curl_2d, compute_divB, average_face_to_cell_B

contains

    subroutine compute_curl_2d(Az, dx, nx, ny, bx, by)
    !
    !   Calculate the discrete curl
    !
    !   Parameters:
    !       - Az : matrix of nodal z-component of magnetic potential
    !       - dx : the cell size
    !       - bx : matrix of cell face x-component magnetic-field
    !       - by : matrix of cell face y-component magnetic-field
    !
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: dx, Az(nx, ny)
        real(8), intent(out) :: bx(nx, ny),  by(nx, ny)
        real(8) :: inv_dx

        ! Calc "bx = dAz/dy" via backward difference
        ! cshift handles periodic BC
        inv_dx = 1.0d0 / dx
        bx =  (Az - cshift(Az, -1, 2)) * inv_dx ! Interior (backward diff. in y)
        by = -(Az - cshift(Az, -1, 1)) * inv_dx ! Interior (backward diff. in x)
        
    end subroutine compute_curl_2d


    subroutine compute_divB(bx, by, dx, nx, ny, divB)
    !
    !   Calculate the discrete divergence of the magnetic field
    !
    !   Parameters:
    !       - dx   : cell size
    !       - bx   : matrix of cell face x-component magnetic-field
    !       - by   : matrix of cell face y-component magnetic-field
    !       - divB : matrix of divergence values for each cell
    !
    !   Notes:
    !       This uses periodic boundary conditions and backward differences:
    !       divB = (bx(i,j) - bx(i-1,j) + by(i,j) - by(i,j-1)) / dx
    !
        integer, intent(in) :: nx, ny
        real(8), intent(in) :: dx, bx(nx, ny), by(nx, ny)
        real(8), intent(out) :: divB(nx, ny)
        real(8) :: inv_dx

        ! Calc divergence in the interior using backward differences
        ! Uses modern Fortran array slicing
        inv_dx = 1.0d0 / dx
        divB(2:nx,2:ny) = ( bx(2:nx,2:ny) - bx(1:nx-1,2:ny) + &   
                            by(2:nx,2:ny) - by(2:nx,1:ny-1) ) * inv_dx

        ! Handle periodic BCs
        ! For i=1: use bx(nx,j) as left neighbor (wrap in x)
        ! For j=1: use by(i,ny) as lower neighbor (wrap in y)
        divB(1,2:ny) = (bx(1,2:ny)-bx(nx,2:ny) +by(1,2:ny)-by(1,1:ny-1))*inv_dx ! Left 
        divB(2:nx,1) = (bx(2:nx,1)-bx(1:nx-1,1)+by(2:nx,1)-by(2:nx,ny))*inv_dx  ! Bottom 
        divB(1,1)    = (bx(1,1) - bx(nx,1) + by(1,1) - by(1,ny))*inv_dx         ! Lower-left

    end subroutine compute_divB


    subroutine average_face_to_cell_B(b_x, b_y, nx, ny, Bx, By)
    !
    !   Calculate the volume-averaged magnetic field
    !       - b_x : matrix of cell face x-component magnetic-field
    !       - b_y : matrix of cell face y-component magnetic-field
    !       - nx  : grid dimensions (x)
    !       - ny  : grid dimensions (y)
    !       - Bx  : matrix of cell Bx (averaged)
    !       - By  : matrix of cell By (averaged)
    !
        integer, intent(in)  :: nx, ny
        real(8), intent(in)  :: b_x(nx, ny), b_y(nx, ny)
        real(8), intent(out) :: Bx(nx, ny),  By(nx, ny)

        ! Avg. x-face fields to cell centers in the x-dir.
        ! For cell i, average b_x at faces i and i-1 (i=1 -> periodic BC)
        Bx(2:nx,:) = 0.5d0 * (b_x(2:nx,:) + b_x(1:nx-1,:)) ! Interior cells
        Bx(1,:)    = 0.5d0 * (b_x(1,:) + b_x(nx,:))        ! Left boundary (wrap)

        ! Average y-face fields to cell centers in the y-dir.
        ! For cell j, average b_y at faces j and j-1 (j=1 -> periodic BC)
        By(:,2:ny) = 0.5d0 * (b_y(:,2:ny) + b_y(:,1:ny-1)) ! Interior cells
        By(:,1)    = 0.5d0 * (b_y(:,1) + b_y(:,ny))        ! Bottom boundary (wrap)

    end subroutine average_face_to_cell_B


end module mhd_field_ops