module mhd_bc
    use mhd_config
    !--------------------------------------------------------------------------
    ! Purpose: Ghost-cell boundary conditions for the 2D MHD solver.
    !
    !   Each call to fill_ghost_cells pads a single N×N physical field into
    !   an (N+2)×(N+2) array.  The mapping is:
    !
    !       physical cell (i, j)  ->  f_pad(i+1, j+1)     i,j ∈ [1,N]
    !       left   ghost          -> f_pad(1,    2:N+1)
    !       right  ghost          ->  f_pad(N+2,  2:N+1)
    !       bottom ghost          ->  f_pad(2:N+1, 1)
    !       top    ghost          ->  f_pad(2:N+1, N+2)
    !
    !   BC types (defined in mhd_config):
    !       BC_PERIODIC - circular wrap (reproduces original cshift behaviour)
    !       BC_OUTFLOW  - zero-gradient: ghost = nearest interior cell
    !       BC_FIXED    - prescribed ambient state via optional f_ambient arg
    !       BC_INFLOW   - driven inflow: same as FIXED; f_ambient must be set
    !                     by the problem setup routine for each relevant field
    !
    !   Usage (main loop, before compute_gradients):
    !       call fill_ghost_cells(rho, rho_pad)
    !       call fill_ghost_cells(P,   P_pad)
    !       ...
    !   For driven inflow on any side:
    !       call fill_ghost_cells(rho, rho_pad, f_ambient=rho_ambient)
    !--------------------------------------------------------------------------

    implicit none
    private
    public fill_ghost_cells

contains

    subroutine fill_ghost_cells(f, f_pad, f_ambient)
    !
    !   Fill an (N+2)x(N+2) padded array from the N×N physical field f,
    !   applying the BC type for each side from mhd_config.
    !
    !   Inputs:
    !       f         - physical field  (N, N)
    !       f_ambient - (optional) prescribed ambient state (N, N);
    !                   required when any side uses BC_FIXED or BC_INFLOW
    !
    !   Output:
    !       f_pad     - padded field  (N+2, N+2)
    !
        real(8), intent(in)  :: f(N, N)
        real(8), intent(out) :: f_pad(N+2, N+2)
        real(8), intent(in), optional :: f_ambient(N, N)

        ! Copy physical domain into interior of padded array 
        f_pad(2:N+1, 2:N+1) = f(1:N, 1:N)

        !--- X-direction ghost columns ---

        ! Left ghost column (xlo)
        select case (bc_xlo)
            case (BC_PERIODIC)
                f_pad(1, 2:N+1) = f(N, 1:N)          ! wrap from right
            case (BC_FIXED, BC_INFLOW)
                if (present(f_ambient)) then
                    f_pad(1, 2:N+1) = f_ambient(1, 1:N)
                else
                    f_pad(1, 2:N+1) = f(1, 1:N)       ! fall back to zero-gradient
                end if
            case default  ! BC_OUTFLOW
                f_pad(1, 2:N+1) = f(1, 1:N)           ! zero-gradient
        end select

        ! Right ghost column (xhi)
        select case (bc_xhi)
            case (BC_PERIODIC)
                f_pad(N+2, 2:N+1) = f(1, 1:N)         ! wrap from left
            case (BC_FIXED, BC_INFLOW)
                if (present(f_ambient)) then
                    f_pad(N+2, 2:N+1) = f_ambient(N, 1:N)
                else
                    f_pad(N+2, 2:N+1) = f(N, 1:N)
                end if
            case default  ! BC_OUTFLOW
                f_pad(N+2, 2:N+1) = f(N, 1:N)
        end select

        !--- Y-direction ghost rows ---

        ! Bottom ghost row (ylo)
        select case (bc_ylo)
            case (BC_PERIODIC)
                f_pad(2:N+1, 1) = f(1:N, N)           ! wrap from top
            case (BC_FIXED, BC_INFLOW)
                if (present(f_ambient)) then
                    f_pad(2:N+1, 1) = f_ambient(1:N, 1)
                else
                    f_pad(2:N+1, 1) = f(1:N, 1)
                end if
            case default  ! BC_OUTFLOW
                f_pad(2:N+1, 1) = f(1:N, 1)
        end select

        ! Top ghost row (yhi)
        select case (bc_yhi)
            case (BC_PERIODIC)
                f_pad(2:N+1, N+2) = f(1:N, 1)         ! wrap from bottom
            case (BC_FIXED, BC_INFLOW)
                if (present(f_ambient)) then
                    f_pad(2:N+1, N+2) = f_ambient(1:N, N)
                else
                    f_pad(2:N+1, N+2) = f(1:N, N)
                end if
            case default  ! BC_OUTFLOW
                f_pad(2:N+1, N+2) = f(1:N, N)
        end select

        ! --- Corners (needed by constrained transport Ez evaluation) ---
        ! Fill last - each corner uses the per-side flags for its two meeting edges.
        ! Fully periodic: diagonal wrap. Mixed or all non-periodic: nearest corner cell.
        call fill_corners(f_pad, f)

    end subroutine fill_ghost_cells


    subroutine fill_corners(f_pad, f)
    !
    !   Fill the four corner ghost cells of f_pad.
    !   Each corner is determined by the BC flags of its two meeting sides.
    !
        real(8), intent(inout) :: f_pad(N+2, N+2)
        real(8), intent(in)    :: f(N, N)

        ! Lower-left  (xlo, ylo)
        f_pad(1,    1)   = corner_val(f, 1, 1, bc_xlo, bc_ylo)
        ! Lower-right (xhi, ylo)
        f_pad(N+2,  1)   = corner_val(f, N, 1, bc_xhi, bc_ylo)
        ! Upper-left  (xlo, yhi)
        f_pad(1,    N+2) = corner_val(f, 1, N, bc_xlo, bc_yhi)
        ! Upper-right (xhi, yhi)
        f_pad(N+2,  N+2) = corner_val(f, N, N, bc_xhi, bc_yhi)

    end subroutine fill_corners


    pure function corner_val(f, ix, iy, bc_this_x, bc_this_y) result(val)
    !
    !   Return the appropriate ghost value for a corner cell.
    !       (ix, iy)     - nearest interior corner indices (1 or N in each dim)
    !       bc_this_x/y  - BC flags for the two sides meeting at this corner
    !
        real(8), intent(in) :: f(N, N)
        integer, intent(in) :: ix, iy, bc_this_x, bc_this_y
        real(8) :: val

        integer :: ix_wrap, iy_wrap

        ! Opposite corner index (1 -> N, N -> 1)
        ix_wrap = N + 1 - ix
        iy_wrap = N + 1 - iy

        if (bc_this_x == BC_PERIODIC .and. bc_this_y == BC_PERIODIC) then
            val = f(ix_wrap, iy_wrap)   ! full diagonal wrap
        else if (bc_this_x == BC_PERIODIC) then
            val = f(ix_wrap, iy)        ! x wraps, y doesn't
        else if (bc_this_y == BC_PERIODIC) then
            val = f(ix, iy_wrap)        ! y wraps, x doesn't
        else
            val = f(ix, iy)             ! nearest interior corner
        end if

    end function corner_val


end module mhd_bc
