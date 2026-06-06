module mhd_field_ops
    use mhd_config

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
    !   Compute the discrete backward-difference curl of Az:
    !       bx(i,j) =  (Az(i,j) - Az(i,j-1)) / dx   = dAz/dy
    !       by(i,j) = -(Az(i,j) - Az(i-1,j)) / dx   = -dAz/dx
    !
    !   For outflow boundaries: interior rows/columns are computed first so
    !   the boundary row/column can copy the nearest interior value (linear
    !   extrapolation of the update). Setting to 0.0 would freeze the
    !   face-centred B at boundary cells and create artificial discontinuities.
    !
        integer, intent(in)  :: nx, ny
        real(8), intent(in)  :: dx, Az(nx, ny)
        real(8), intent(out) :: bx(nx, ny), by(nx, ny)
        real(8) :: inv_dx

        inv_dx = 1.0d0 / dx

        ! --- bx = dAz/dy (backward diff in y) ---
        ! Interior first so j=2 is ready for the ylo boundary copy
        bx(:, 2:ny) = (Az(:, 2:ny) - Az(:, 1:ny-1)) * inv_dx

        ! ylo boundary (j=1)
        select case (bc_ylo)
            case (BC_PERIODIC)
                bx(:, 1) = (Az(:, 1) - Az(:, ny)) * inv_dx
            case default  ! outflow: copy nearest interior update (linear extrapolation)
                bx(:, 1) = bx(:, 2)
        end select

        ! --- by = -dAz/dx (backward diff in x) ---
        ! Interior first so i=2 is ready for the xlo boundary copy
        by(2:nx, :) = -(Az(2:nx, :) - Az(1:nx-1, :)) * inv_dx

        ! xlo boundary (i=1)
        select case (bc_xlo)
            case (BC_PERIODIC)
                by(1, :) = -(Az(1, :) - Az(nx, :)) * inv_dx
            case default  ! outflow: copy nearest interior update
                by(1, :) = by(2, :)
        end select

    end subroutine compute_curl_2d


    subroutine compute_divB(bx, by, dx, nx, ny, divB)
    !
    !   Compute the discrete backward-difference divergence of B:
    !       divB(i,j) = (bx(i,j) - bx(i-1,j) + by(i,j) - by(i,j-1)) / dx
    !
    !   Boundary cells use BC flags from mhd_config.
    !   For outflow: zero-gradient ghost -> b_ghost = b_boundary -> diff = 0.
    !   Note: divB is a diagnostic only; it does not feed back into the solver.
    !
        integer, intent(in)  :: nx, ny
        real(8), intent(in)  :: dx, bx(nx, ny), by(nx, ny)
        real(8), intent(out) :: divB(nx, ny)
        real(8) :: inv_dx

        inv_dx = 1.0d0 / dx

        ! Interior cells
        divB(2:nx, 2:ny) = (bx(2:nx,2:ny) - bx(1:nx-1,2:ny) + &
                            by(2:nx,2:ny) - by(2:nx,1:ny-1)) * inv_dx

        ! Left column (i=1), interior rows
        select case (bc_xlo)
            case (BC_PERIODIC)
                divB(1,2:ny) = (bx(1,2:ny) - bx(nx,2:ny) + &
                                by(1,2:ny) - by(1,1:ny-1)) * inv_dx
            case default  ! outflow: bx(0,j) = bx(1,j) -> x-contribution = 0
                divB(1,2:ny) = (by(1,2:ny) - by(1,1:ny-1)) * inv_dx
        end select

        ! Bottom row (j=1), interior columns
        select case (bc_ylo)
            case (BC_PERIODIC)
                divB(2:nx,1) = (bx(2:nx,1) - bx(1:nx-1,1) + &
                                by(2:nx,1) - by(2:nx,ny)) * inv_dx
            case default  ! outflow: by(i,0) = by(i,1) -> y-contribution = 0
                divB(2:nx,1) = (bx(2:nx,1) - bx(1:nx-1,1)) * inv_dx
        end select

        ! Corner (i=1, j=1)
        select case (bc_xlo)
            case (BC_PERIODIC)
                ! x-contribution wraps
                select case (bc_ylo)
                    case (BC_PERIODIC)
                        divB(1,1) = (bx(1,1)-bx(nx,1) + by(1,1)-by(1,ny)) * inv_dx
                    case default
                        divB(1,1) = (bx(1,1)-bx(nx,1)) * inv_dx
                end select
            case default
                select case (bc_ylo)
                    case (BC_PERIODIC)
                        divB(1,1) = (by(1,1)-by(1,ny)) * inv_dx
                    case default
                        divB(1,1) = 0.0d0
                end select
        end select

    end subroutine compute_divB


    subroutine average_face_to_cell_B(b_x, b_y, nx, ny, Bx, By)
    !
    !   Average face-centred B fields to cell centres.
    !   For cell i: Bx(i) = 0.5*(b_x(i) + b_x(i-1))
    !   For cell j: By(j) = 0.5*(b_y(j) + b_y(j-1))
    !
    !   Boundary cells use BC flags from mhd_config.
    !   For outflow: zero-gradient ghost -> b_ghost = b_boundary,
    !   so Bx(1) = b_x(1) and By(:,1) = b_y(:,1).
    !
        integer, intent(in)  :: nx, ny
        real(8), intent(in)  :: b_x(nx, ny), b_y(nx, ny)
        real(8), intent(out) :: Bx(nx, ny),  By(nx, ny)

        ! Interior x-averages
        Bx(2:nx, :) = 0.5d0 * (b_x(2:nx,:) + b_x(1:nx-1,:))

        ! xlo boundary (i=1)
        select case (bc_xlo)
            case (BC_PERIODIC)
                Bx(1,:) = 0.5d0 * (b_x(1,:) + b_x(nx,:))
            case default  ! outflow: ghost = cell -> average collapses to cell value
                Bx(1,:) = b_x(1,:)
        end select

        ! Interior y-averages
        By(:, 2:ny) = 0.5d0 * (b_y(:,2:ny) + b_y(:,1:ny-1))

        ! ylo boundary (j=1)
        select case (bc_ylo)
            case (BC_PERIODIC)
                By(:,1) = 0.5d0 * (b_y(:,1) + b_y(:,ny))
            case default
                By(:,1) = b_y(:,1)
        end select

    end subroutine average_face_to_cell_B


end module mhd_field_ops