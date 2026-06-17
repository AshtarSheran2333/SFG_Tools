submodule (FRAME_READERS) GRO_READER
    use, intrinsic :: iso_fortran_env, only: output_unit
    use PRETTY_PRINT
    use UTILS_ERROR
    implicit none

#include "utils_error_macros.h"

    character(*), parameter :: gro_format = "(i5,2a5,i5,3f8.3,3f8.3)"

    contains
    
    module procedure gro_read_header
        implicit none
        logical :: is_open
        integer :: ierr
        
        res = .false.

        inquire(this%file, opened = is_open)

        if(.not. is_open) return
        
        read(this%file, "(A)", iostat = ierr) !comment line
        if(ierr .ne. 0) return

        read(this%file, *, iostat = ierr) this%frame%n_atoms !number of atoms
        if(ierr .ne. 0) return

        res = .true.
    end procedure gro_read_header

    module procedure gro_open_file
        implicit none
        integer :: ierr
        logical :: is_open
        integer(int32) :: res_index, atom_index
        character(len = 5) :: res_name, atom_name
        real(kind=real64), dimension(3) :: dummy

        res = -1

        if(present(filename1)) error_stop("Calling gro reader with two arguments is not allowed")
        
        if(allocated(this%frame%positions)) then
            deallocate(this%frame%positions)
        end if
        if(allocated(this%frame%velocities)) then
            deallocate(this%frame%velocities)
        end if
        if(allocated(this%frame%names)) then
            deallocate(this%frame%names)
        end if
        
        this%frame%frame_number = 0
        this%frame%n_atoms = 0
        this%frame%has_velocities = .false.
        
        inquire(this%file, opened = is_open)
        if(is_open) then
            close(this%file)
        end if
        
        write(output_unit,'( "Opening ", A, " file...")') trim(filename)
        write(output_unit,f_line) ""

        open(newunit = this%file, file = filename, status = 'old', iostat = ierr)
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        
        if(.not. this%read_header()) return
        
        this%prev_n_atoms = this%frame%n_atoms

        allocate(this%frame%positions(3,this%frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        allocate(this%frame%velocities(3,this%frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        allocate(this%frame%names(this%frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return

        !attempt to read velocities
        this%frame%has_velocities = .true.
        read(this%file, gro_format, iostat = ierr) &
            res_index, res_name, atom_name, atom_index, dummy, dummy
        if(ierr .ne. 0) then
            this%frame%has_velocities = .false.
            !attempt to read without velocities
            read(this%file, gro_format, iostat = ierr) &
                res_index, res_name, atom_name, atom_index, dummy
            if(ierr .ne. 0) then
                !problem
                res = ierr
                return
            end if
        end if
        
        rewind(this%file, iostat = ierr) !rewind
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
        inquire(this%file, opened = is_open)
        if(.not. is_open) return

        res = -2; if(.not. this%read_header()) return
        
        !check wheather the number of atoms changed, if yes, reallocate this%atoms
        if(this%frame%n_atoms .ne. this%prev_n_atoms) then
            error_stop("gro reader does not support variable number of atoms") 
            !here one can implement variable number of atoms, but I dont want this functionality
        end if
        
        if(this%frame%has_velocities) then
            read(this%file, gro_format, iostat = ierr) &
                (res_index, res_name, this%frame%names(i), atom_index, this%frame%positions(:,i), this%frame%velocities(:,i), &
                i = 1, this%frame%n_atoms)
        else
            read(this%file, gro_format, iostat = ierr) &
                (res_index, res_name, this%frame%names(i), atom_index, this%frame%positions(:,i), &
                i = 1, this%frame%n_atoms)
        end if
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
            
        !convert the positions
        this%frame%positions = this%frame%positions * nm_to_angstrom
        !convert the velocities
        if(this%frame%has_velocities) this%frame%velocities = this%frame%velocities * nm_ps_to_a_fs

        !read the box size - just discard it
        read(this%file, *, iostat = ierr) box
        res = -3; if(ierr .ne. 0) return
        
        this%frame%frame_number = this%frame%frame_number + 1
        res = 0;
    end procedure gro_read_frame

    module procedure gro_skip_frame
        implicit none
        integer :: ierr
        logical :: is_open
        character(len=128) :: dummy
        integer(int64) :: i

        res = -1
        inquire(this%file, opened = is_open)
        if(.not. is_open) return

        res = -2; if(.not. this%read_header()) return
        
        read(this%file, "(A)", iostat = ierr) &
            (dummy , i = 1, this%frame%n_atoms+1) !empty read +1 for the box
        if(ierr .ne. 0) then
            res = ierr
            return
        end if

        this%frame%frame_number = this%frame%frame_number + 1
        res = 0;
    end procedure gro_skip_frame

    module procedure gro_rewind_file
        implicit none
        logical :: is_open
        integer :: ierr
        
        inquire(this%file, opened = is_open)
        res = -1; if(.not. is_open) return

        rewind(this%file, iostat = ierr)
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        
        this%frame%frame_number = 0
        res = 0
    end procedure gro_rewind_file

    module procedure gro_close_file
        implicit none
        logical :: is_open
        
        inquire(this%file, opened = is_open)
        res = -1; if(.not. is_open) return
        
        close(this%file, iostat = res)
        if(res .ne. 0) return

        this%file = 0

        if(allocated(this%frame%positions)) deallocate(this%frame%positions)
        if(allocated(this%frame%velocities)) deallocate(this%frame%velocities)
        if(allocated(this%frame%names)) deallocate(this%frame%names)
        
        this%frame%n_atoms = 0
        this%frame%frame_number = 0
        this%frame%has_velocities = .false.
        res = 0
    end procedure gro_close_file

    module procedure gro_is_open
        implicit none
        inquire(this%file, opened = res)
    end procedure gro_is_open

end submodule GRO_READER
