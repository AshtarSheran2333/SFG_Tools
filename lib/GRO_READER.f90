submodule (FRAME_READERS) GRO_READER
    use, intrinsic :: iso_fortran_env, only: output_unit
    use UTILS_ERROR
    implicit none

#include "utils_error_macros.h"

    character(*), parameter :: gro_format = "(i5,2a5,i5,3f8.3,3f8.3)"

    contains
    
    function read_header() result(res)
        implicit none
        logical :: res
        logical :: is_open
        integer :: ierr
        
        res = .false.

        inquire(fr_file, opened = is_open)

        if(.not. is_open) return
        
        read(fr_file, "(A)", iostat = ierr) !comment line
        if(ierr .ne. 0) return

        read(fr_file, *, iostat = ierr) fr_frame%n_atoms !number of atoms
        if(ierr .ne. 0) return

        res = .true.
    end function read_header

    module procedure gro_open_file
        implicit none
        integer :: ierr
        logical :: is_open
        integer(int32) :: res_index, atom_index
        character(len = 5) :: res_name, atom_name
        real(kind=real64), dimension(3) :: dummy

        res = -1

        if(present(filename1)) error_stop("Calling gro reader with two arguments is not allowed")
        
        if(allocated(fr_frame%positions)) then
            deallocate(fr_frame%positions)
        end if
        if(allocated(fr_frame%velocities)) then
            deallocate(fr_frame%velocities)
        end if
        if(allocated(fr_frame%names)) then
            deallocate(fr_frame%names)
        end if
        
        fr_frame%frame_number = 0
        fr_frame%n_atoms = 0
        fr_frame%has_velocities = .false.
        
        inquire(fr_file, opened = is_open)
        if(is_open) then
            close(fr_file)
        end if
        
        write(output_unit,'( "Opening ", A, " file...")') trim(filename)

        open(newunit = fr_file, file = filename, status = 'old', iostat = ierr)
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        
        if(.not. read_header()) return
        
        prev_n_atoms = fr_frame%n_atoms

        allocate(fr_frame%positions(3,fr_frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        allocate(fr_frame%velocities(3,fr_frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        allocate(fr_frame%names(fr_frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return

        !attempt to read velocities
        fr_frame%has_velocities = .true.
        read(fr_file, gro_format, iostat = ierr) &
            res_index, res_name, atom_name, atom_index, dummy, dummy
        if(ierr .ne. 0) then
            fr_frame%has_velocities = .false.
            !attempt to read without velocities
            read(fr_file, gro_format, iostat = ierr) &
                res_index, res_name, atom_name, atom_index, dummy
            if(ierr .ne. 0) then
                !problem
                res = ierr
                return
            end if
        end if
        
        rewind(fr_file, iostat = ierr) !rewind
        if(ierr .ne. 0) return

        res = 0
    end procedure gro_open_file

    module procedure gro_read_frame
        implicit none
        integer(int64) :: i
        integer(int32) :: res_index, atom_index
        real(real32), dimension(3) :: box
        character(len = 5) :: res_name, atom_name
        integer :: ierr
        logical :: is_open
        
        res = -1
        inquire(fr_file, opened = is_open)
        if(.not. is_open) return

        res = -2; if(.not. read_header()) return
        
        !check wheather the number of atoms changed, if yes, reallocate fr_atoms
        if(fr_frame%n_atoms .ne. prev_n_atoms) then
            error_stop("gro reader does not support variable number of atoms") 
            !here one can implement variable number of atoms, but I dont want this functionality
        end if
        
        if(fr_frame%has_velocities) then
            read(fr_file, gro_format, iostat = ierr) &
                (res_index, res_name, fr_frame%names(i), atom_index, fr_frame%positions(:,i), fr_frame%velocities(:,i), &
                i = 1, fr_frame%n_atoms)
        else
            read(fr_file, gro_format, iostat = ierr) &
                (res_index, res_name, fr_frame%names(i), atom_index, fr_frame%positions(:,i), &
                i = 1, fr_frame%n_atoms)
        end if
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
            
        !convert the positions
        fr_frame%positions = fr_frame%positions * nm_to_angstrom
        !convert the velocities
        if(fr_frame%has_velocities) fr_frame%velocities = fr_frame%velocities * nmpps_to_hartree

        !read the box size - just discard it
        read(fr_file, *, iostat = ierr) box
        res = -3; if(ierr .ne. 0) return
        
        fr_frame%frame_number = fr_frame%frame_number + 1
        res = 0;
    end procedure gro_read_frame

    module procedure gro_skip_frame
        implicit none
        integer :: ierr
        logical :: is_open
        character(len=128) :: dummy
        integer(int64) :: i

        res = -1
        inquire(fr_file, opened = is_open)
        if(.not. is_open) return

        res = -2; if(.not. read_header()) return
        
        read(fr_file, "(A)", iostat = ierr) &
            (dummy , i = 1, fr_frame%n_atoms+1) !empty read +1 for the box
        if(ierr .ne. 0) then
            res = ierr
            return
        end if

        fr_frame%frame_number = fr_frame%frame_number + 1
        res = 0;
    end procedure gro_skip_frame

    module procedure gro_rewind_file
        implicit none
        logical :: is_open
        integer :: ierr
        
        inquire(fr_file, opened = is_open)
        res = -1; if(.not. is_open) return

        rewind(fr_file, iostat = ierr)
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        
        fr_frame%frame_number = 0
        res = 0
    end procedure gro_rewind_file

    module procedure gro_close_file
        implicit none
        logical :: is_open
        
        inquire(fr_file, opened = is_open)
        res = -1; if(.not. is_open) return
        
        close(fr_file, iostat = res)
        if(res .ne. 0) return

        fr_file = 0

        if(allocated(fr_frame%positions)) deallocate(fr_frame%positions)
        if(allocated(fr_frame%velocities)) deallocate(fr_frame%velocities)
        if(allocated(fr_frame%names)) deallocate(fr_frame%names)
        
        fr_frame%n_atoms = 0
        fr_frame%frame_number = 0
        fr_frame%has_velocities = .false.
        res = 0
    end procedure gro_close_file

    module procedure gro_is_open
        implicit none
        inquire(fr_file, opened = res)
    end procedure gro_is_open

end submodule GRO_READER
