program SFG_BINDER
use iso_fortran_env
use BINDER_UI
use FRAME_READERS
use BOXDATA
use INSTANTANEOUS_SURFACE
use SFG_STRUCTURE
!TODO make a binder module - responsible for holding the binder, writing, reading...

implicit none

#include "utils_error_macros.h"

type(boxdata_type) ::                                       bd

class(frame_reader), allocatable ::                         fr

type(instantaneous_surface_type) ::                         instasurf

type(sfg_structure_type) ::                                 struct

integer(int64) ::											step

!##########################################EVALUATE PROGRAM SWITCHES###########################################

call evaluate_program_options(fr)

!TODO read BOXDATA

!TODO READ INTERFACE

!TODO READ STRUCT

!TODO INIT BINDER

!TODO MAIN LOOP

do step = 1, bd%NSTEP

	!TODO reading frame
	
	!TODO estimate center of the water slab
	
    !TODO read interface each INTERFACE_SKIP
	
	!TODO loop over all chromophores in the system, assign it a layer
    
    !TODO write binder frame
    
    !TODO print progress
end do

!TODO finalize() cleanup

!TODO print done
print*, ""

contains

end program SFG_BINDER
