program SFG_BINDER
use iso_fortran_env
use BINDER_UI
use FRAME_READERS
use BOXDATA
use INSTANTANEOUS_SURFACE
use BINDER_FILE
use SFG_STRUCTURE
!TODO make a binder module - responsible for holding the binder, writing, reading...

implicit none

#include "utils_error_macros.h"

type(boxdata_type) ::                                       bd

class(frame_reader), allocatable ::                         fr

type(instantaneous_surface_type) ::                         instasurf

type(sfg_structure_type) ::                                 struct

type(binder_type) ::                                        binder

integer(int64) ::                                           step

integer ::  res

real(real64) ::                                             z_center

!##########################################EVALUATE PROGRAM SWITCHES###########################################

call evaluate_program_options(fr)

call bd%read_boxdata()

!init instasurf & its files
error_io_check(instasurf%open_bin_file(read_only = .true., must_exist = .true., name = "interface.bin"), "Problem opening interface.bin file")

!read structure
error_io_check(struct%read_structure("struct.txt"), "unable to read structure file")

!TODO INIT BINDER
error_io_check(binder%init(struct), "unable to init binder")
error_io_check(binder%open_file(), "unable to open binder file")

write(output_unit,f_line) heading(wavy_pattern, "Binder calculation started")
write(output_unit,f_line) ""

!################################ MAIN LOOP ####################################
do step = 1, bd%NSTEP

    !read frame
    error_io_check(fr%read_frame(), "unable to read a frame")
    
    if( (mod(step, bd%INTERFACE_SKIP) == 1) .or. (bd%INTERFACE_SKIP == 1) ) then
        !TODO better reporting - need to report step number
        error_io_check(instasurf%read_next(), "unable to read instantaneous surface")
        !estimate center of the liquid slab - geometric average of the instantaneous surfaces
        z_center = instasurf%get_center()
    end if
    
    !The binder is initialized by the structure, obrain the whole binder for all the groups specified by the structure
    error_io_check(binder%fill_binder(fr, struct, instasurf, bd, z_center), "binder assignment error")
    
    !write binder frame
    error_io_check(binder%write_frame(), "unable to write binder frame")

    call print_main_loop_progress(step, bd%NSTEP)

end do !end of the main loop

!cleanup
call finalize()

write(output_unit,f_line) ""
write(output_unit,f_line) heading(wavy_pattern, "Binder calculation - DONE")

contains

subroutine finalize()
    implicit none
    
    call instasurf%close_bin_file()
    error_io_check(binder%close_file(), "unable to close binder frame")
    
end subroutine finalize

end program SFG_BINDER
