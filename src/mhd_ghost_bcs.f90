module mhd_ghost_BCs
    use mhd_config, only: N
    implicit none
    private
    public fill_outflow_ghosts 

contains

    subroutine fill_outflow_ghosts(f_in, f_ghost_out, ghost_width, direction)
        real(8), intent(in)  :: f_in(N,N)
        real(8), intent(out) :: f_ghost_out(:,:)
        integer, intent(in) :: ghost_width
        character(len=*), intent(in) :: direction
        
        integer :: i,j,g, gx,gy
        real(8) :: damp
        
        gx = N + 2*ghost_width
        gy = N + 2*ghost_width
        
        ! Copy physical domain
        f_ghost_out(ghost_width+1:ghost_width+N, ghost_width+1:ghost_width+N) = f_in
        
        ! Fill ghost zones
        if (direction=='x' .or. direction=='both') then
            ! Left ghosts
            do g=1,ghost_width
                damp = 1.0d0 - 0.6d0*real(g,8)
                do j=1,N
                    f_ghost_out(g,j) = f_in(1,j) * damp
                end do
            end do
            ! Right ghosts  
            do g=1,ghost_width
                damp = 1.0d0 - 0.6d0*real(g,8)
                do j=1,N
                    i = N + ghost_width + g
                    f_ghost_out(i,j) = f_in(N,j) * damp
                end do
            end do
        end if
        
        if (direction=='y' .or. direction=='both') then
            ! Bottom ghosts
            do g=1,ghost_width
                damp = 1.0d0 - 0.6d0*real(g,8)
                do i=1,N
                    f_ghost_out(i,g) = f_in(i,1) * damp
                end do
            end do
            ! Top ghosts
            do g=1,ghost_width
                damp = 1.0d0 - 0.6d0*real(g,8)
                do i=1,N
                    j = N + ghost_width + g
                    f_ghost_out(i,j) = f_in(i,N) * damp
                end do
            end do
        end if
        
    end subroutine fill_outflow_ghosts


end module mhd_ghost_BCs