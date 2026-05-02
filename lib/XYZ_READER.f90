submodule (FRAME_READERS) XYZ_READER
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    use UTILS_ERROR
    implicit none
    
#include "utils_error_macros.h"

    contains
    
    function read_header_with_velocities() result(res)
        implicit none
        logical :: res
        logical :: is_open
        integer :: ierr
        integer(int32) :: vel_n_atoms, pos_n_atoms
        character(len = 256) :: line

        res = .false.
        
        inquire(fr_file, opened = is_open)
        if(.not. is_open) return
        inquire(fr_file1, opened = is_open)
        if(.not. is_open) return

        read(fr_file, *, iostat = ierr) pos_n_atoms
        if(ierr .ne. 0) return

        read(fr_file1, *, iostat = ierr) vel_n_atoms
        if(ierr .ne. 0) return

        if(pos_n_atoms .ne. vel_n_atoms) then
            write(error_unit, "(A)") "Error number of atoms does not match for xyz positions and velocities file"
            return
        end if
        
        fr_frame%n_atoms = pos_n_atoms

        read(fr_file, "(A)", iostat = ierr) !comment line
        read(fr_file1, "(A)", iostat = ierr) !comment line

        if(ierr .ne. 0) return
        
        res = .true.
    end function read_header_with_velocities
    
    function read_header_no_velocities() result(res)
        implicit none
        logical :: res
        logical :: is_open
        integer :: ierr

        res = .false.

        inquire(fr_file, opened = is_open)
        if(.not. is_open) return

        read(fr_file, *, iostat = ierr) fr_frame%n_atoms !number of atoms
        if(ierr .ne. 0) return

        read(fr_file, "(A)", iostat = ierr) !comment line
        if(ierr .ne. 0) return

        res = .true. 
    end function read_header_no_velocities

    function open_with_velocities(posfile, velfile) result(res)
        implicit none
        character(*), intent(in) :: posfile, velfile
        integer :: res
        integer :: ierr
        integer(int32) :: vel_n_atoms, pos_n_atoms

        res = -1
        
        open(newunit = fr_file, file = posfile, status = 'old', iostat = ierr)
        if(ierr .ne. 0) then
            write(error_unit, "(A)") "ERROR: opening the positions xyz file"
            res = ierr
            return
        end if

        open(newunit = fr_file1, file = velfile, status = 'old', iostat = ierr)
        if(error_unit .ne. 0) then
            write(output_unit, "(A)") "ERROR: opening the velocities xyz file"
            res = ierr
            return
        end if
        
        read(fr_file, *, iostat = ierr) pos_n_atoms
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        read(fr_file1, *, iostat = ierr) vel_n_atoms
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        
        if(pos_n_atoms .ne. vel_n_atoms) then
            write(error_unit, "(A)") "ERROR: number of atoms does not match for xyz positions and velocities file"
            return
        end if
        
        fr_frame%n_atoms = pos_n_atoms
        prev_n_atoms = fr_frame%n_atoms

        rewind(fr_file, iostat = ierr) !rewind
        if(ierr .ne. 0) return
        rewind(fr_file1, iostat = ierr) !rewind
        if(ierr .ne. 0) return

        res = 0
    end function open_with_velocities
    
    function open_no_velocities(posfile) result(res)
        implicit none
        character(*), intent(in) :: posfile
        integer :: res
        integer :: ierr

        res = -1
        
        open(newunit = fr_file, file = posfile, status = 'old', iostat = ierr)
        if(ierr .ne. 0) then
            write(error_unit, "(A)") "ERROR: opening the positions xyz file"
            res = ierr
            return
        end if
        
        read(fr_file, *, iostat = ierr) fr_frame%n_atoms
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        
        prev_n_atoms = fr_frame%n_atoms

        rewind(fr_file, iostat = ierr) !rewind
        if(ierr .ne. 0) return

        res = 0
    end function open_no_velocities

    module procedure xyz_open_file
        implicit none
        integer :: ierr
        logical :: is_open
        integer(int32) :: res_index, atom_index
        character(len = 5) :: res_name, atom_name

        res = -1

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

        inquire(fr_file1, opened = is_open)
        if(is_open) then
            close(fr_file)
        end if

        write(output_unit,'( "Opening ", A, " file...")') trim(filename)

        if(present(filename1)) then
            res = open_with_velocities(filename, filename1)
            if(res .eq. 0) then
                fr_frame%has_velocities = .true.
            end if
        else
            res = open_no_velocities(filename)
        end if
        
        allocate(fr_frame%positions(3,fr_frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        allocate(fr_frame%velocities(3,fr_frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        allocate(fr_frame%names(fr_frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return

        res = 0
    end procedure xyz_open_file

    module procedure xyz_read_frame
        implicit none
        integer(int64) :: i
        integer :: ierr

        res = -1
        if(fr_frame%has_velocities) then
            if(.not. read_header_with_velocities()) return 
        else
            if(.not. read_header_no_velocities()) return
        end if

        !check wheather the number of atoms changed, if yes, reallocate fr_atoms
        if(fr_frame%n_atoms .ne. prev_n_atoms) then
            error_stop("xyz reader does not support variable number of atoms") 
            !here one can implement variable number of atoms, but I dont want this functionality
        end if

        !TODO - now expecting that the atomnames match... (it should be checked)
        if(fr_frame%has_velocities) then
            read(fr_file, *, iostat = ierr) (fr_frame%names(i), fr_frame%positions(:,i), i = 1, fr_frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
            read(fr_file1, *, iostat = ierr) (fr_frame%names(i), fr_frame%velocities(:,i), i = 1, fr_frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
        else
            read(fr_file, *, iostat = ierr) (fr_frame%names(i), fr_frame%positions(:,i), i = 1, fr_frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
        end if

        fr_frame%frame_number = fr_frame%frame_number + 1
        res = 0
        return
    end procedure xyz_read_frame

    module procedure xyz_skip_frame
        implicit none
        integer(int64) :: i
        integer :: ierr
        character(len=128) :: dummy

        res = -1
        if(fr_frame%has_velocities) then
            if(.not. read_header_with_velocities()) return 
        else
            if(.not. read_header_no_velocities()) return
        end if

        !check wheather the number of atoms changed, if yes, reallocate fr_atoms
        if(fr_frame%n_atoms .ne. prev_n_atoms) then
            error_stop("xyz reader does not support variable number of atoms") 
            !here one can implement variable number of atoms, but I dont want this functionality
        end if

        !TODO - now expecting that the atomnames match... (it should be checked)
        if(fr_frame%has_velocities) then
            read(fr_file, "(A)", iostat = ierr) (dummy, i = 1, fr_frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
            read(fr_file1, "(A)", iostat = ierr) (dummy, i = 1, fr_frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
        else
            read(fr_file, "(A)", iostat = ierr) (dummy, i = 1, fr_frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
        end if

        fr_frame%frame_number = fr_frame%frame_number + 1
        res = 0
    end procedure xyz_skip_frame

    module procedure xyz_rewind_file
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

        inquire(fr_file1, opened = is_open)
        if(is_open) then
            rewind(fr_file1, iostat = ierr)
            if(ierr .ne. 0) then
                res = ierr
                return
            end if
        end if
        
        fr_frame%frame_number = 0
        res = 0
    end procedure xyz_rewind_file

    module procedure xyz_close_file
        implicit none
        logical :: is_open
        integer :: ierr
        
        inquire(fr_file, opened = is_open)
        res = -1; if(.not. is_open) return

        close(fr_file, iostat = ierr)
        if(ierr .ne. 0) return
        fr_file = 0

        inquire(fr_file1, opened = is_open)
        if(is_open) then
            close(fr_file1, iostat = ierr)
            if(ierr .ne. 0) return
            fr_file1 = 0
        end if

        if(allocated(fr_frame%positions)) deallocate(fr_frame%positions)
        if(allocated(fr_frame%velocities)) deallocate(fr_frame%velocities)
        if(allocated(fr_frame%names)) deallocate(fr_frame%names)
        
        fr_frame%n_atoms = 0
        fr_frame%frame_number = 0
        fr_frame%has_velocities = .false.
        res = 0 
    end procedure xyz_close_file

    module procedure xyz_is_open
        implicit none
        logical :: is_open
        
        res = .false.
        
        inquire(fr_file, opened = is_open)
        if(.not. is_open) return 
        
        if(fr_frame%has_velocities) then 
            inquire(fr_file1, opened = is_open)
            if(.not. is_open) return
        end if

        res = .true.
    end procedure xyz_is_open

end submodule XYZ_READER
