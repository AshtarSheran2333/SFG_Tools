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

!TODO cleanup of the garbage variables - most of them
type(boxdata_type) ::                                       bd

class(frame_reader), pointer ::                             fr
type(trr_frame_reader), allocatable, target ::              trr
type(gro_frame_reader), allocatable, target ::              gro
type(xyz_frame_reader), allocatable, target ::              xyz

real*8, dimension (:,:,:), allocatable ::                   points
!points replaces lx,ly,lz --- coordinates of both interfaces

real*8, dimension(:,:), allocatable ::                      density_function
!holds up all the density functions (distance,density,updensity,downdensity)
real*8, dimension(:), allocatable ::						density_axis

real*8, dimension(3) ::                                     prev_rod_pos,&
                                                            rod_pos,&
                                                            rod_start,&
                                                            diff
!rod_pos replaces l1, l2, l3 --- coordinates of "poking rod"
!diff replaces xdiff, ydiff, zdiff --- stores temporary distance of oxygen atoms from "poking rod"
!rod_start replaces fx, fy, fz --- coordinates of rod starting position

integer*8, dimension (:,:,:), allocatable ::                prev_index

integer*8, dimension(:,:), allocatable ::                   ij_array

integer*8, dimension(2,2) ::                                ij_est

integer*8, dimension(3) ::                                  n_points
!n_points replaces nlx, nly, nlz --- number of divisions in each direction given by dividing boxdimensions by dl

integer, dimension(2) ::                                    uji,&
                                                            dji

integer ::                                                  ierr,&
                                                            f_interface,&
                                                            f_grid_interface,&
                                                            f_interfacebin

real*8 ::                                                   r,&
                                                            p,&
                                                            pdiff,&
                                                            u,&
                                                            urmin,&
                                                            drmin,&
                                                            uzdiff,&
                                                            dzdiff,&
                                                            zcenter,&
                                                            search_radius,&
                                                            prev_pdiff
!r - distance of given oxygen atom from the end of "poking rod"
!p - sum[oxygen atoms] if condition is met - exp(-r**2/(2*E**2))/((2*pi*E**2)**1.5)
!u
!urmin - min of !up! r = ||(xdiff,ydiff,zdiff)||
!drmin - min of !down! r = ||(xdiff,ydiff,zdiff)||
!uzdiff 
!dzdiff
!zcenter - z coordinate approximately in the middle of water slab

integer*8 ::                                                step,&
                                                            d,&
                                                            i,&
                                                            j,&
                                                            k,&
                                                            m,&
                                                            s,&
                                                            nop,&
                                                            temp
!step - iterator through steps
!d - sumcheck of points each frame -> should be allways n_points(1)*n_points(2)*2
!i,j,k -general iterators
!m -iterates through NO
!s -temporary storage for writing step number to the output file

logical ::                                                  found_uper_interface,&
                                                            found_bottom_interface,&
                                                            vmdout = .FALSE.,&
                                                            vmdpbc = .TRUE.,&
                                                            bin_header_done = .FALSE.,&
                                                            nocheck = .false.
!found_uper_interface/min holds up info if given point met conditions to become rod_max/min  

character(len = 256) ::                                     errmsg

character(len=60), allocatable ::                           arg,&
                                                            arg1

real*8, parameter ::                                        &
                                                            E=2.4,&
                                                            tollerance=0.004,&
                                                            waterDensity=0.03336 !N/A^3

type(instantaneous_surface_type) :: instasurf

type(density_profile_type) :: heavy_atoms_density
type(density_profile_type) :: waters_density
type(density_profile_type), dimension(:), allocatable :: up_group_densities, bot_group_densities

real(real64), dimension(2) :: is_ret
!##########################################EVALUATE PROGRAM SWITCHES###########################################

call evaluate_program_options()

!TODO deal with frame reader
select case(ui_filetype)
case (ui_filetype_gro)
    allocate(gro)
    fr => gro
    error_io_check(fr%open_file(ui_filename1), "unable to open file")
case (ui_filetype_trr)
    allocate(trr)
    fr => trr
    error_io_check(fr%open_file(ui_filename1, ui_filename2), "unable to open file")
case default
    !TODO other filetypes
    error_stop("unknown input filetype")
end select


call bd%read_boxdata()

! TODO init instasurf & its files
call instasurf%init(bd)

error_io_check(instasurf%write_grid_interface(), "unable to write grid")
error_io_check(instasurf%open_xyz_file(), "unable to open interface.xyz")
error_io_check(instasurf%open_bin_file(), "unable to open interface.bin")

error_io_check(read_struct("struct.txt"), "unable to read structure file")
!TODO init density
call heavy_atoms_density%init(bd%DENSITY_R_START, bd%DENSITY_R_END, bd%DENSITY_N_POINTS, bd%BOX_DIMENSIONS, "heavy_atoms")
call waters_density%init(bd%DENSITY_R_START, bd%DENSITY_R_END, bd%DENSITY_N_POINTS, bd%BOX_DIMENSIONS, "waters")

call init_group_densities()

!TODO main loop

do step = 1, bd%NSTEP, bd%INTERFACE_SKIP
    
    print "(a,a,i,$)", char(13), "STEP: ", step
    
    error_io_check(fr%read_frame(), "unable to read a frame") !reading frame

    !TODO calculate interface each INTERFACE_SKIP frames
    error_io_check(instasurf%calculate(fr%frame, bd), "cannot calculate instasurf")
    !TODO or read instantaneous surface from a file
    
    !TODO calculation of densities & water dipole moments...

    !density of liquid heavy atoms
    associate(dp => heavy_atoms_density%bins(:))
    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(bd, fr, instasurf, heavy_atoms_density) &
    !$omp private(i, is_ret) &
    !$omp reduction(+:dp)
    do i = 1, size(bd%LIQUID_HEAVY_ATOMS)
        
        is_ret = instasurf%get_distances(fr%frame%positions(:,i), bd)
        call heavy_atoms_density%add_point(is_ret(1), 1.0_real64)

    end do
    !$omp end parallel do
    end associate
    call heavy_atoms_density%next_frame()

    !density of waters
    associate(dp => waters_density%bins(:))
    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(bd, fr, instasurf, waters_density, sfg_structure) &
    !$omp private(i, j, diff, is_ret) &
    !$omp reduction(+:dp)
    do i = 1, size(sfg_structure(1)%sfg_units)

        !get "center of mass" (geometric average of all bases)
        diff = 0
        do j = 1, size(sfg_structure(1)%sfg_units(i)%chromophores)
            diff = diff + pbc_wrap(fr%frame%positions(:,sfg_structure(1)%sfg_units(i)%chromophores(1)%base), bd)
        end do
        diff = diff /  size(sfg_structure(1)%sfg_units(i)%chromophores)

        is_ret = instasurf%get_distances(diff, bd)
        call waters_density%add_point(is_ret(1), 1.0_real64)

    end do
    !$omp end parallel do
    end associate
    call waters_density%next_frame()

    do j = 1, size(sfg_structure)
        call get_group_density(j)
    end do

    !TODO writing the instantaneous surface frames
    error_io_check(instasurf%write_xyz_frame(), "did not write the xyz frame")
    error_io_check(instasurf%write_bin_frame(), "did not write the bin frame")

end do

!write the density files
call heavy_atoms_density%write_file()
call waters_density%write_file()

do j = 1, size(sfg_structure)
    call up_group_densities(j)%write_file()
    call bot_group_densities(j)%write_file()
end do

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

    !density of structgroup
    associate(udp => up_group_densities(index)%bins(:), &
                bdp => bot_group_densities(index)%bins(:))
    !$omp parallel do &
    !$omp default(none) &
    !$omp shared(bd, fr, instasurf, index, up_group_densities, bot_group_densities, sfg_structure) &
    !$omp private(i, j, diff, is_ret) &
    !$omp reduction(+:udp, bdp)
    do i = 1, size(sfg_structure(index)%sfg_units)

        !get "center of mass" (geometric average of all bases)
        diff = 0
        do j = 1, sfg_structure(index)%sfg_units(i)%n_unique_bases
            diff = diff + pbc_wrap(fr%frame%positions(:,sfg_structure(index)%sfg_units(i)%unique_bases(j)), bd)
        end do
        diff = diff /  sfg_structure(index)%sfg_units(i)%n_unique_bases

        is_ret = instasurf%get_distances(diff, bd)
        call up_group_densities(index)%add_point(is_ret(2), 1.0_real64)
        call bot_group_densities(index)%add_point(is_ret(1), 1.0_real64)

    end do
    !$omp end parallel do
    end associate
    call up_group_densities(index)%next_frame()
    call bot_group_densities(index)%next_frame()
end subroutine

end program SFG_INTERFACE
