module mhd_write_h5
    !--------------------------------------------------------------------
    ! Purpose: Write the simulated primitive / conserved field data to 
    !          a single HDF5 output file with physics parameters
    !--------------------------------------------------------------------
    use hdf5
    implicit none
    private
    public write_prims_to_hdf5

contains

    subroutine write_prims_to_hdf5(filename, rho_dat, P_dat, Bx_dat, By_dat, Vx_dat, Vy_dat, &
                                   time_dat, gamma_val, M_s_val, beta_val)
        implicit none

        character(len=*), intent(in) :: filename
        character(len=5) :: int_string
        real(kind=8), dimension(:,:,:), intent(in) :: rho_dat, P_dat, Bx_dat, By_dat, Vx_dat, Vy_dat
        real(kind=8), dimension(:), intent(in) :: time_dat
        real(kind=8), intent(in) :: gamma_val, M_s_val, beta_val
        
        integer :: error
        integer(HSIZE_T) :: data_dims(2), time_dims(1), scalar_dims(1)
        integer(HID_T)   :: file_id, group_id1, group_id2, group_id3, group_id4, group_id5
        integer(HID_T)   :: group_id6, group_id7, group_id_params
        integer(HID_T)   :: dspace_id, dset_id
        integer :: i

        ! Initialize error
        error = 0  

        ! Open interface
        call h5open_f(error)

        ! Open file (create / overwrite)
        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
        if (error /= 0) then
            print *, "Error opening file"
            stop
        end if

        ! Create groups for each field
        call h5gcreate_f(file_id, "rho", group_id1, error)
        call h5gcreate_f(file_id, "P",  group_id2, error)
        call h5gcreate_f(file_id, "Bx", group_id3, error)
        call h5gcreate_f(file_id, "By", group_id4, error)
        call h5gcreate_f(file_id, "Vx", group_id5, error)
        call h5gcreate_f(file_id, "Vy", group_id6, error)
        call h5gcreate_f(file_id, "time", group_id7, error)
        call h5gcreate_f(file_id, "parameters", group_id_params, error)
        if (error /= 0) then
            print *, "Error creating groups"
            stop
        end if


        ! Loop over snapshots (existing code unchanged)
        do i = 1, size(rho_dat, 1)

            data_dims = [size(rho_dat, 2), size(rho_dat, 3)]
            call h5screate_simple_f(2, data_dims, dspace_id, error)
            if (error /= 0) then
                print *, "Error creating dataspace"
                stop
            end if

            write(int_string, '(I0)') i  

            !--- Create dataset for rho ---
            call h5dcreate_f(group_id1, "SNAPSHOT"//trim(int_string), H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
            if (error /= 0) then
                print *, "Error creating dataset for rho"
                stop
            end if
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, rho_dat(i, :, :), data_dims, error)
            if (error /= 0) then
                print *, "Error writing data for rho"
                stop
            end if
            call h5dclose_f(dset_id, error)

            !--- Create dataset for P ---
            call h5dcreate_f(group_id2, "SNAPSHOT"//trim(int_string), H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
            if (error /= 0) then
                print *, "Error creating dataset for P"
                stop
            end if
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, P_dat(i, :, :), data_dims, error)
            if (error /= 0) then
                print *, "Error writing data for P"
                stop
            end if
            call h5dclose_f(dset_id, error)

            !--- Create dataset for Bx ---
            call h5dcreate_f(group_id3, "SNAPSHOT"//trim(int_string), H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
            if (error /= 0) then
                print *, "Error creating dataset for Bx"
                stop
            end if
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Bx_dat(i, :, :), data_dims, error)
            if (error /= 0) then
                print *, "Error writing data for Bx"
                stop
            end if
            call h5dclose_f(dset_id, error)

            !--- Create dataset for By ---
            call h5dcreate_f(group_id4, "SNAPSHOT"//trim(int_string), H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
            if (error /= 0) then
                print *, "Error creating dataset for By"
                stop
            end if
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, By_dat(i, :, :), data_dims, error)
            if (error /= 0) then
                print *, "Error writing data for By"
                stop
            end if
            call h5dclose_f(dset_id, error)

            !--- Create dataset for Vx ---
            call h5dcreate_f(group_id5, "SNAPSHOT"//trim(int_string), H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
            if (error /= 0) then
                print *, "Error creating dataset for Vx"
                stop
            end if
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Vx_dat(i, :, :), data_dims, error)
            if (error /= 0) then
                print *, "Error writing data for Vx"
                stop
            end if
            call h5dclose_f(dset_id, error)

            !--- Create dataset for Vy ---
            call h5dcreate_f(group_id6, "SNAPSHOT"//trim(int_string), H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
            if (error /= 0) then
                print *, "Error creating dataset for Vy"
                stop
            end if
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, Vy_dat(i, :, :), data_dims, error)
            if (error /= 0) then
                print *, "Error writing data for Vy"
                stop
            end if
            call h5dclose_f(dset_id, error)

            call h5sclose_f(dspace_id, error)

        end do

        !--- Write 1D time array ---
        time_dims = [size(time_dat)]
        call h5screate_simple_f(1, time_dims, dspace_id, error)
        call h5dcreate_f(group_id7, "sim_time", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time_dat, time_dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dspace_id, error)

        ! Write physics parameters as scalar datasets
        scalar_dims(1) = 1
        call h5screate_simple_f(1, scalar_dims, dspace_id, error)
        
        ! Write gamma
        call h5dcreate_f(group_id_params, "gamma", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, gamma_val, scalar_dims, error)
        call h5dclose_f(dset_id, error)
        
        ! Write M_s
        call h5dcreate_f(group_id_params, "M_s", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, M_s_val, scalar_dims, error)
        call h5dclose_f(dset_id, error)
        
        ! Write beta
        call h5dcreate_f(group_id_params, "beta", H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, beta_val, scalar_dims, error)
        call h5dclose_f(dset_id, error)
        
        call h5sclose_f(dspace_id, error)

        ! Close groups and file
        call h5gclose_f(group_id1, error)
        call h5gclose_f(group_id2, error)
        call h5gclose_f(group_id3, error)
        call h5gclose_f(group_id4, error)
        call h5gclose_f(group_id5, error)
        call h5gclose_f(group_id6, error)
        call h5gclose_f(group_id7, error)
        call h5gclose_f(group_id_params, error)
        call h5fclose_f(file_id, error)

        call h5close_f(error)

    end subroutine write_prims_to_hdf5

end module mhd_write_h5