program SFG_INTERFACE

use iso_fortran_env
use INTERFACE_UI
use FRAME_READERS
use BOXDATA
use INSTANTANEOUS_SURFACE
use DENSITY
use SFG_STRUCTURE

implicit none

#include "utils_error_macros.h"

type(boxdata_type) ::                                       bd

class(frame_reader), allocatable ::                         fr

type(instantaneous_surface_type) ::                         instasurf

type(sfg_structure_type) ::                                 struct

type(density_profile_type), dimension(:), allocatable ::    up_group_densities,&
                                                            bot_group_densities

integer(int64) ::                                           step,&
                                                            sk
integer ::                                                  gr

!############################ BEGIN INITIALIZATION #############################

!evaluate program options and open trajectory...
call evaluate_program_options(fr)

call bd%read_boxdata()

!init instasurf & its files
call instasurf%init(bd)
error_io_check(instasurf%write_grid_interface(), "unable to write grid")

if(ui_vmd_out) then
    error_io_check(instasurf%open_xyz_file(), "unable to open interface.xyz")
end if

!TODO implement skip of the interface calculation - read from a file
error_io_check(instasurf%open_bin_file(), "unable to open interface.bin")

!read the struct file, so we can evaluate densities of the groups...
error_io_check(struct%read_structure("struct.txt"), "unable to read structure file")

!init densities
call init_group_densities()

write(output_unit,f_line) heading(wavy_pattern, "Interface calculation started")
write(output_unit,f_line) ""

!################################ MAIN LOOP ####################################
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
    do gr = 1, size(struct%groups)
        call get_group_density(gr)
    end do

    call print_main_loop_progress(step, bd%NSTEP)

    !reads the skipped frames, evaluates densities of the skipped frames
    do sk = 1, min(bd%INTERFACE_SKIP-1, bd%NSTEP-step)
        error_io_check(fr%read_frame(), "unable to read a frame") !reading frame
        
        do gr = 1, size(struct%groups)
            call get_group_density(gr)
        end do

        call print_main_loop_progress(step+sk, bd%NSTEP)
    end do

end do !end of the main loop

!finalize the analyses
call finalize()

write(output_unit,f_line) ""
write(output_unit,f_line) heading(wavy_pattern, "Interface calculation - DONE")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS AND SUBROUTINES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

!TODO those subroutines touch global variables... it should not be problem, but I would be careful about that!

subroutine init_group_densities()
    implicit none
    integer :: i
    
    if(allocated(up_group_densities)) deallocate(up_group_densities)
    allocate(up_group_densities(size(struct%groups)))
    if(allocated(bot_group_densities)) deallocate(bot_group_densities)
    allocate(bot_group_densities(size(struct%groups)))

    do i = 1, size(struct%groups)
        call up_group_densities(i)%init(bd%DENSITY_R_START, bd%DENSITY_R_END, bd%DENSITY_N_POINTS, bd%BOX_DIMENSIONS, "up_"//trim(adjustl(struct%groups(i)%name)))
        call bot_group_densities(i)%init(bd%DENSITY_R_START, bd%DENSITY_R_END, bd%DENSITY_N_POINTS, bd%BOX_DIMENSIONS, "bot_"//trim(adjustl(struct%groups(i)%name)))
    end do
end subroutine init_group_densities

subroutine get_group_density(index)
    implicit none
    integer, intent(in) :: index
    integer :: i, j
    real(real64), dimension(2) :: is_ret
    real(real64), dimension(3) :: com

    !TODO the openmp approach like this is permitted by ifx, but highly nonstandard
    !density of structgroup
    !associate(udp => up_group_densities(index)%bins(:), &
    !            bdp => bot_group_densities(index)%bins(:))
    !!$omp parallel do &
    !!$omp default(none) &
    !!$omp shared(bd, fr, instasurf, index, up_group_densities, bot_group_densities, struct) &
    !!$omp private(i, j, com, is_ret) &
    !!$omp reduction(+:udp, bdp)
    do i = 1, size(struct%groups(index)%sfg_units)

        !get "center of mass" (geometric average of all bases)
        com = 0
        do j = 1, struct%groups(index)%sfg_units(i)%n_unique_bases
            com = com + pbc_wrap(fr%frame%positions(:,struct%groups(index)%sfg_units(i)%unique_bases(j)), bd)
        end do
        com = com / struct%groups(index)%sfg_units(i)%n_unique_bases

        is_ret = instasurf%get_distances(com, bd)
        call up_group_densities(index)%add_point(is_ret(2), 1.0_real64)
        call bot_group_densities(index)%add_point(is_ret(1), 1.0_real64)

    end do
    !!$omp end parallel do
    !end associate
    call up_group_densities(index)%next_frame()
    call bot_group_densities(index)%next_frame()
end subroutine

subroutine finalize()
    implicit none
    integer :: gr

    if(ui_vmd_out) then
        call instasurf%close_bin_file()
    end if

    call instasurf%close_bin_file()
    
    do gr = 1, size(struct%groups)
        call up_group_densities(gr)%write_file()
        call bot_group_densities(gr)%write_file()
    end do
end subroutine finalize

end program SFG_INTERFACE
