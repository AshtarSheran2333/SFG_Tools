submodule (FRAME_READERS) XYZ_READER
    use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
    use PRETTY_PRINT
    use UTILS_ERROR
    implicit none
    
#include "utils_error_macros.h"

    contains
    
    module procedure xyz_read_header_with_velocities
        implicit none
        logical :: is_open
        integer :: ierr
        integer(int32) :: vel_n_atoms, pos_n_atoms
        character(len = 256) :: line

        res = .false.
        
        inquire(this%file, opened = is_open)
        if(.not. is_open) return
        inquire(this%file1, opened = is_open)
        if(.not. is_open) return

        read(this%file, *, iostat = ierr) pos_n_atoms
        if(ierr .ne. 0) return

        read(this%file1, *, iostat = ierr) vel_n_atoms
        if(ierr .ne. 0) return

        if(pos_n_atoms .ne. vel_n_atoms) then
            write(error_unit, "(A)") "Error number of atoms does not match for xyz positions and velocities file"
            return
        end if
        
        this%frame%n_atoms = pos_n_atoms

        read(this%file, "(A)", iostat = ierr) !comment line
        read(this%file1, "(A)", iostat = ierr) !comment line

        if(ierr .ne. 0) return
        
        res = .true.
    end procedure xyz_read_header_with_velocities
    
    module procedure xyz_read_header_no_velocities
        implicit none
        logical :: is_open
        integer :: ierr

        res = .false.

        inquire(this%file, opened = is_open)
        if(.not. is_open) return

        read(this%file, *, iostat = ierr) this%frame%n_atoms !number of atoms
        if(ierr .ne. 0) return

        read(this%file, "(A)", iostat = ierr) !comment line
        if(ierr .ne. 0) return

        res = .true. 
    end procedure xyz_read_header_no_velocities

    module procedure xyz_open_with_velocities
        implicit none
        integer :: ierr
        integer(int32) :: vel_n_atoms, pos_n_atoms

        res = -1
        
        open(newunit = this%file, file = posfile, status = 'old', iostat = ierr)
        if(ierr .ne. 0) then
            write(error_unit, "(A)") "ERROR: opening the positions xyz file"
            res = ierr
            return
        end if

        open(newunit = this%file1, file = velfile, status = 'old', iostat = ierr)
        if(error_unit .ne. 0) then
            write(output_unit, "(A)") "ERROR: opening the velocities xyz file"
            res = ierr
            return
        end if
        
        read(this%file, *, iostat = ierr) pos_n_atoms
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        read(this%file1, *, iostat = ierr) vel_n_atoms
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        
        if(pos_n_atoms .ne. vel_n_atoms) then
            write(error_unit, "(A)") "ERROR: number of atoms does not match for xyz positions and velocities file"
            return
        end if
        
        this%frame%n_atoms = pos_n_atoms
        this%prev_n_atoms = this%frame%n_atoms

        rewind(this%file, iostat = ierr) !rewind
        if(ierr .ne. 0) return
        rewind(this%file1, iostat = ierr) !rewind
        if(ierr .ne. 0) return

        res = 0
    end procedure xyz_open_with_velocities
    
    module procedure xyz_open_no_velocities
        implicit none
        integer :: ierr

        res = -1
        
        open(newunit = this%file, file = posfile, status = 'old', iostat = ierr)
        if(ierr .ne. 0) then
            write(error_unit, "(A)") "ERROR: opening the positions xyz file"
            res = ierr
            return
        end if
        
        read(this%file, *, iostat = ierr) this%frame%n_atoms
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        
        this%prev_n_atoms = this%frame%n_atoms

        rewind(this%file, iostat = ierr) !rewind
        if(ierr .ne. 0) return

        res = 0
    end procedure xyz_open_no_velocities

    module procedure xyz_open_file
        implicit none
        integer :: ierr
        logical :: is_open
        integer(int32) :: res_index, atom_index
        character(len = 5) :: res_name, atom_name

        res = -1

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

        inquire(this%file1, opened = is_open)
        if(is_open) then
            close(this%file)
        end if

        write(output_unit,'( "Opening ", A, " file...")') trim(filename)
        write(output_unit,f_line) ""

        if(present(filename1)) then
            res = this%open_with_velocities(filename, filename1)
            if(res .eq. 0) then
                this%frame%has_velocities = .true.
            end if
        else
            res = this%open_no_velocities(filename)
        end if
        
        allocate(this%frame%positions(3,this%frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        allocate(this%frame%velocities(3,this%frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        allocate(this%frame%names(this%frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return

        res = 0
    end procedure xyz_open_file

    module procedure xyz_read_frame
        implicit none
        integer(int64) :: i
        integer :: ierr

        res = -1
        if(this%frame%has_velocities) then
            if(.not. this%read_header_with_velocities()) return 
        else
            if(.not. this%read_header_no_velocities()) return
        end if

        !check wheather the number of atoms changed, if yes, reallocate this%atoms
        if(this%frame%n_atoms .ne. this%prev_n_atoms) then
            error_stop("xyz reader does not support variable number of atoms") 
            !here one can implement variable number of atoms, but I dont want this functionality
        end if

        !TODO - now expecting that the atomnames match... (it should be checked)
        !TODO the implicit loops might not work properly...
        if(this%frame%has_velocities) then
            read(this%file, *, iostat = ierr) (this%frame%names(i), this%frame%positions(:,i), i = 1, this%frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
            read(this%file1, *, iostat = ierr) (this%frame%names(i), this%frame%velocities(:,i), i = 1, this%frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
        else
            read(this%file, *, iostat = ierr) (this%frame%names(i), this%frame%positions(:,i), i = 1, this%frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
        end if

        this%frame%velocities = this%frame%velocities * hartree_velocity_to_a_fs
        this%frame%frame_number = this%frame%frame_number + 1
        res = 0
        return
    end procedure xyz_read_frame

    module procedure xyz_skip_frame
        implicit none
        integer(int64) :: i
        integer :: ierr
        character(len=128) :: dummy

        res = -1
        if(this%frame%has_velocities) then
            if(.not. this%read_header_with_velocities()) return 
        else
            if(.not. this%read_header_no_velocities()) return
        end if

        !check wheather the number of atoms changed, if yes, reallocate this%atoms
        if(this%frame%n_atoms .ne. this%prev_n_atoms) then
            error_stop("xyz reader does not support variable number of atoms") 
            !here one can implement variable number of atoms, but I dont want this functionality
        end if

        !TODO - now expecting that the atomnames match... (it should be checked)
        if(this%frame%has_velocities) then
            read(this%file, "(A)", iostat = ierr) (dummy, i = 1, this%frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
            read(this%file1, "(A)", iostat = ierr) (dummy, i = 1, this%frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
        else
            read(this%file, "(A)", iostat = ierr) (dummy, i = 1, this%frame%n_atoms)
            res = ierr; if(ierr .ne. 0) return
        end if

        this%frame%frame_number = this%frame%frame_number + 1
        res = 0
    end procedure xyz_skip_frame

    module procedure xyz_rewind_file
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

        inquire(this%file1, opened = is_open)
        if(is_open) then
            rewind(this%file1, iostat = ierr)
            if(ierr .ne. 0) then
                res = ierr
                return
            end if
        end if
        
        this%frame%frame_number = 0
        res = 0
    end procedure xyz_rewind_file

    module procedure xyz_close_file
        implicit none
        logical :: is_open
        integer :: ierr
        
        inquire(this%file, opened = is_open)
        res = -1; if(.not. is_open) return

        close(this%file, iostat = ierr)
        if(ierr .ne. 0) return
        this%file = 0

        inquire(this%file1, opened = is_open)
        if(is_open) then
            close(this%file1, iostat = ierr)
            if(ierr .ne. 0) return
            this%file1 = 0
        end if

        if(allocated(this%frame%positions)) deallocate(this%frame%positions)
        if(allocated(this%frame%velocities)) deallocate(this%frame%velocities)
        if(allocated(this%frame%names)) deallocate(this%frame%names)
        
        this%frame%n_atoms = 0
        this%frame%frame_number = 0
        this%frame%has_velocities = .false.
        res = 0 
    end procedure xyz_close_file

    module procedure xyz_is_open
        implicit none
        logical :: is_open
        
        res = .false.
        
        inquire(this%file, opened = is_open)
        if(.not. is_open) return 
        
        if(this%frame%has_velocities) then 
            inquire(this%file1, opened = is_open)
            if(.not. is_open) return
        end if

        res = .true.
    end procedure xyz_is_open

end submodule XYZ_READER
