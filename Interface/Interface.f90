program SFG_INTERFACE

use iso_fortran_env
use INTERFACE_UI
use FRAME_READERS
use BOXDATA
use INSTANTANEOUS_SURFACE
use DENSITY
use SFG_STRUCT

implicit none

#include "utils_error_macros.h"

type(boxdata_type) ::                                       bd

class(frame_reader), allocatable ::                         fr

integer(int64) ::                                           step,&
                                                            sk,&
                                                            gr

type(instantaneous_surface_type) ::                         instasurf

type(density_profile_type), dimension(:), allocatable ::    up_group_densities,&
                                                            bot_group_densities

!##########################################EVALUATE PROGRAM SWITCHES###########################################

call evaluate_program_options()

!allocate the correct frame reader todo can be moved to UI_CORE or something like that...
select case(ui_filetype)
case (ui_filetype_gro)
    allocate(gro_frame_reader :: fr)
    error_io_check(fr%open_file(ui_filename1), "unable to open "//ui_filename1)
case (ui_filetype_trr)
    allocate(trr_frame_reader :: fr)
    error_io_check(fr%open_file(ui_filename1, ui_filename2), "unable to open files "//ui_filename1//" "//ui_filename2)
case (ui_filetype_xyz)
    allocate(xyz_frame_reader :: fr)
    error_io_check(fr%open_file(ui_filename1, ui_filename2), "unable to open files "//ui_filename1//" "//ui_filename2)
case default
    error_stop("unknown input filetype")
end select

call bd%read_boxdata()

!init instasurf & its files
call instasurf%init(bd)
error_io_check(instasurf%write_grid_interface(), "unable to write grid")

!TODO open xyz only when vmdout
error_io_check(instasurf%open_xyz_file(), "unable to open interface.xyz")
!TODO implement skip of the interface calculation - read from a file
error_io_check(instasurf%open_bin_file(), "unable to open interface.bin")

!read the struct file, so we can evaluate densities of the groups...
error_io_check(read_struct("struct.txt"), "unable to read structure file")

!init densities
call init_group_densities()

write(output_unit,f_line) heading(wavy_pattern, "Interface calculation started")
write(output_unit,f_line) ""

!main loop
do step = 1, bd%NSTEP, bd%INTERFACE_SKIP

    error_io_check(fr%read_frame(), "unable to read a frame") !reading frame

    !calculate the instantaneous surface each INTERFACE_SKIP frames
    error_io_check(instasurf%calculate(fr%frame, bd), "cannot calculate instasurf")
    !TODO or read instantaneous surface from a file

    !writing the instantaneous surface frames
    error_io_check(instasurf%write_bin_frame(), "did not write the bin frame")
    if(ui_vmd_out) then
        error_io_check(instasurf%write_xyz_frame(), "did not write the xyz frame")
    end if
    
    !calculation of densities
    do gr = 1, size(sfg_structure)
        call get_group_density(gr)
    end do

    !reads the skipped frames, evaluates densities of the skipped frames
    do sk = 1, min(bd%INTERFACE_SKIP-1, step-bd%INTERFACE_SKIP-1)
        error_io_check(fr%read_frame(), "unable to read a frame") !reading frame
        
        do gr = 1, size(sfg_structure)
            call get_group_density(gr)
        end do
    end do

    call print_main_loop_progress(step, bd%NSTEP)

end do !end of the main loop

!finalize the analyses
call finalize()

write(output_unit,f_line) ""
write(output_unit,f_line) heading(wavy_pattern, "Interface calculation - DONE")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS AND SUBROUTINES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

subroutine init_group_densities()
    implicit none
    !TODO wrap SFG_struct into a type....
    integer :: i
    
    if(allocated(up_group_densities)) deallocate(up_group_densities)
    allocate(up_group_densities(size(sfg_structure)))
    if(allocated(bot_group_densities)) deallocate(bot_group_densities)
    allocate(bot_group_densities(size(sfg_structure)))

    do i = 1, size(sfg_structure)
        call up_group_densities(i)%init(bd%DENSITY_R_START, bd%DENSITY_R_END, bd%DENSITY_N_POINTS, bd%BOX_DIMENSIONS, "up_"//trim(adjustl(sfg_structure(i)%name)))
        call bot_group_densities(i)%init(bd%DENSITY_R_START, bd%DENSITY_R_END, bd%DENSITY_N_POINTS, bd%BOX_DIMENSIONS, "bot_"//trim(adjustl(sfg_structure(i)%name)))
    end do
end subroutine init_group_densities

subroutine get_group_density(index)
    implicit none
    !TODO wrap SFG_struct into a type, pass it as an argument
    integer, intent(in) :: index
    integer :: i, j
    real(real64), dimension(2) :: is_ret
    real(real64), dimension(3) :: com

    !density of structgroup
    associate(udp => up_group_densities(index)%bins(:), &
                bdp => bot_group_densities(index)%bins(:))
    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(bd, fr, instasurf, index, up_group_densities, bot_group_densities, sfg_structure) &
    !$omp private(i, j, com, is_ret) &
    !$omp reduction(+:udp, bdp)
    do i = 1, size(sfg_structure(index)%sfg_units)

        !get "center of mass" (geometric average of all bases)
        com = 0
        do j = 1, sfg_structure(index)%sfg_units(i)%n_unique_bases
            com = com + pbc_wrap(fr%frame%positions(:,sfg_structure(index)%sfg_units(i)%unique_bases(j)), bd)
        end do
        com = com / sfg_structure(index)%sfg_units(i)%n_unique_bases

        is_ret = instasurf%get_distances(com, bd)
        call up_group_densities(index)%add_point(is_ret(2), 1.0_real64)
        call bot_group_densities(index)%add_point(is_ret(1), 1.0_real64)

    end do
    !$omp end parallel do
    end associate
    call up_group_densities(index)%next_frame()
    call bot_group_densities(index)%next_frame()
end subroutine

subroutine finalize()
    implicit none
    integer :: gr

    do gr = 1, size(sfg_structure)
        call up_group_densities(gr)%write_file()
        call bot_group_densities(gr)%write_file()
    end do
end subroutine finalize

end program SFG_INTERFACE
