submodule (FRAME_READERS) TRR_READER
    use, intrinsic :: iso_fortran_env, only: output_unit
    use UTILS_ERROR
    implicit none

#include "utils_error_macros.h"
    
    contains
    
    module procedure trr_read_header
        implicit none
        integer(int32) :: skip = 0
        character(:), allocatable :: cw
        logical :: is_open
        integer :: ierr
        
        res = .false.

        inquire(this%file, opened = is_open)

        if(.not. is_open) return
        
        !magic number
        read(this%file, iostat = ierr) skip !1993
        if(ierr .ne. 0) return

        if(skip /= 1993) then
            write(output_unit,*) "corrupted trr file header"
            return
        end if

        read(this%file, iostat = ierr) skip !max str len
        if(ierr .ne. 0) return

        read(this%file, iostat = ierr) skip !actual str len
        if(ierr .ne. 0) return

        allocate(character(len=skip) :: cw)
        
        read(this%file, iostat = ierr) cw !"GMX_trn_file"
        if(ierr .ne. 0) return

        if(cw .ne. "GMX_trn_file") then
            write(output_unit,*) "unexpected trr string"
            return
        end if

        deallocate(cw)
        
        !àead the whole header in one go
        associate(h => this%header)
        read(this%file, iostat = ierr) &
            h%sizes%ir_size,    &
            h%sizes%e_size,     &
            h%sizes%box_size,   &
            h%sizes%vir_size,   &
            h%sizes%pres_size,  &
            h%sizes%top_size,   &
            h%sizes%sym_size,   &
            h%sizes%pos_size,   &
            h%sizes%vel_size,   &
            h%sizes%forces_size,&
            h%n_atoms,          &
            h%step_number,      &
            h%nre,              &
            h%sim_time,         &
            h%lambda 
        end associate
        if(ierr .ne. 0) return

        res = .true.
    end procedure trr_read_header
    
    module procedure trr_read_ir !not implemented - just skipping the data
        implicit none
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(this%header%sizes%ir_size .eq. 0) then 
            res = .true.
            return !nothing to read
        end if

        inquire(this%file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(this%file, pos = position + this%header%sizes%ir_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong
        
        res = .true.
    end procedure trr_read_ir

    module procedure trr_read_e !not implemented - just skipping the data
        implicit none
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(this%header%sizes%e_size .eq. 0) then 
            res = .true.
            return !nothing to read
        end if
        
        inquire(this%file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(this%file, pos = position + this%header%sizes%e_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong
        
        res = .true.
    end procedure trr_read_e

    module procedure trr_read_box !not implemented - just skipping the data
        implicit none
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(this%header%sizes%box_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(this%file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(this%file, pos = position + this%header%sizes%box_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end procedure trr_read_box

    module procedure trr_read_vir !not implemented - just skipping the data
        implicit none
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(this%header%sizes%vir_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(this%file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(this%file, pos = position + this%header%sizes%vir_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end procedure trr_read_vir

    module procedure trr_read_pres !not implemented - just skipping the data
        implicit none
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(this%header%sizes%pres_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(this%file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(this%file, pos = position + this%header%sizes%pres_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end procedure trr_read_pres

    module procedure trr_read_top !not implemented - just skipping the data
        implicit none
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(this%header%sizes%top_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(this%file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(this%file, pos = position + this%header%sizes%top_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end procedure trr_read_top
    
    module procedure trr_read_sym !not implemented - just skipping the data
        implicit none
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(this%header%sizes%sym_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(this%file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(this%file, pos = position + this%header%sizes%sym_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end procedure trr_read_sym
    
    module procedure trr_read_positions
        implicit none
        integer :: ierr
        real(real32), dimension(3) :: read_sp
        integer(int32) :: element_sz
        integer(int32) :: i
        logical :: is_open
        
        res = .false.
        
        if(this%header%sizes%pos_size .eq. 0) then
            return !nothing to read
        end if

        inquire(this%file, opened = is_open)
        if(.not. is_open) return
        
        element_sz = this%header%sizes%pos_size / (3 * this%header%n_atoms)
        
        if( (element_sz * 3 * this%header%n_atoms) .ne. this%header%sizes%pos_size ) &
            error_stop("position size is not congruent with the number of atoms")
        
        if( this%header%n_atoms .ne. size(this%frame%positions,2) ) &
            error_stop("trr reader internal arrays size does not match the trr header n_atoms")
            
        if( (element_sz .ne. 4) .and. (element_sz .ne. 8) ) &
            error_stop("trr positions are neither single nor double precision")
        
        if(element_sz .eq. 8) then
            !read double precision
            read(this%file, iostat = ierr) (this%frame%positions(:,i), i = 1, this%header%n_atoms) 
        else
            !read single precision, store in double precision
            do i = 1, this%header%n_atoms
                  read(this%file, iostat = ierr) read_sp
                  if(ierr .ne. 0) exit
                  this%frame%positions(:,i) = read_sp
            end do
        end if

        if(ierr .ne. 0) return !something went wrong

        !convert the positions
        this%frame%positions = this%frame%positions * nm_to_angstrom
        
        res = .true.
    end procedure trr_read_positions
    
    module procedure trr_read_velocities
        implicit none
        integer :: ierr
        real(real32), dimension(3) :: read_sp
        integer(int32) :: element_sz
        integer(int32) :: i
        logical :: is_open
        
        res = .false.
        
        if(this%header%sizes%vel_size .eq. 0) then
            if(this%frame%has_velocities) return !error, config changed
        end if

        inquire(this%file, opened = is_open)
        if(.not. is_open) return
        
        element_sz = this%header%sizes%vel_size / (3 * this%header%n_atoms)
        
        if( (element_sz * 3 * this%header%n_atoms) .ne. this%header%sizes%vel_size ) &
            error_stop("trr velocity size is not congruent with the number of atoms")
        
        if( this%header%n_atoms .ne. size(this%frame%velocities,2) ) &
            error_stop("trr reader internal n_atoms does not match the trr header n_atoms")
        
        if( (element_sz .ne. 4) .and. (element_sz .ne. 8) ) &
            error_stop("trr velocities are neither single nor double precision")
        
        if(element_sz .eq. 8) then
            !read double precision
            read(this%file, iostat = ierr) (this%frame%velocities(:,i), i = 1, this%header%n_atoms) 
        else
            !read single precision, store in double precision
            do i = 1, this%header%n_atoms
                  read(this%file, iostat = ierr) read_sp
                  if(ierr .ne. 0) exit
                  this%frame%velocities(:,i) = read_sp
            end do
        end if

        if(ierr .ne. 0) return !something went wrong

        !convert the velocities
        this%frame%velocities = this%frame%velocities * nmpps_to_hartree

        this%frame%has_velocities = .true.
        res = .true.
    end procedure trr_read_velocities

    module procedure trr_read_forces !not implemented - just skipping the data
        implicit none
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(this%header%sizes%forces_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(this%file, pos = position, opened = is_open)
        
        if(.not. is_open) return
        
        read(this%file, pos = position + this%header%sizes%forces_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end procedure trr_read_forces

    module procedure trr_fill_atom_names
        implicit none
        integer :: ierr
        integer(int64) :: natoms
        character(128) :: dummy_str
        real(real32) :: dummy_real
        integer(int32) dummy_int        
        
        res = .false.

        write(output_unit,'( "Opening ", A, " file...")') trim(filename)

        open(newunit = this%file1, file = filename, status = 'old', iostat = ierr)
        if(ierr .ne. 0) return

        read(this%file1, "(A)", iostat = ierr) dummy_str !read header str
        if(ierr .ne. 0) return

        read(this%file1, *, iostat = ierr) natoms !read natoms
        if(ierr .ne. 0) return
        
        if(natoms .ne. this%frame%n_atoms) return
        
        if(.not. allocated(this%frame%names)) return !this%atoms must be allocated
        if(size(this%frame%names) .ne. this%frame%n_atoms) return !this%atoms size must correspond to the file
        
        !fill in the atoms
        read(this%file1, "(I5,2A5)", iostat = ierr) &
            (dummy_int, dummy_str, this%frame%names(natoms),&
            natoms = 1, this%frame%n_atoms)

        if(ierr .ne. 0) return
        
        close(this%file1)
        
        res = .true.
    end procedure trr_fill_atom_names
    
    module procedure trr_open_file
        implicit none
        integer :: ierr, file1
        logical :: is_open
        res = -1

        ! not necessary, but I want to have the atomnames!
        ! - should work even if the following line is commented out
        if(.not. present(filename1)) error_stop("trr_open_file must be called with both arguments")
        
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

        open(newunit = this%file, file = filename, form = "unformatted",&
                access = 'stream', convert = 'big_endian', status = 'old', iostat = ierr)
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        
        if(.not. this%read_header()) return
        
        this%frame%n_atoms = this%header%n_atoms
        this%prev_n_atoms = this%frame%n_atoms

        if(this%header%sizes%vel_size > 0) this%frame%has_velocities = .true.

        allocate(this%frame%positions(3,this%frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        allocate(this%frame%velocities(3,this%frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        allocate(this%frame%names(this%frame%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        
        rewind(this%file, iostat = ierr) !rewind
        if(ierr .ne. 0) return

        if(.not. this%fill_atom_names(filename1)) return

        res = 0
    end procedure trr_open_file

    module procedure trr_read_frame
        implicit none
        logical :: is_open
        
        res = -1
        inquire(this%file, opened = is_open)
        if(.not. is_open) return

        res = -2; if(.not. this%read_header()) return
        
        !check wheather the number of atoms changed, if yes, reallocate this%atoms
        if(this%header%n_atoms .ne. this%prev_n_atoms) then
            error_stop("trr reader does not support variable number of atoms") 
            !here one can implement variable number of atoms, but I dont want this functionality
        end if
        
        res = -3; if(.not. this%read_ir()) return
        res = -4; if(.not. this%read_e()) return
        res = -5; if(.not. this%read_box()) return
        res = -6; if(.not. this%read_vir()) return
        res = -7; if(.not. this%read_pres()) return
        res = -8; if(.not. this%read_top()) return
        res = -9; if(.not. this%read_sym()) return
        res = -10; if(.not. this%read_positions()) return
        res = -11; if(.not. this%read_velocities()) return
        res = -12; if(.not. this%read_forces()) return

        this%frame%frame_number = this%frame%frame_number + 1
        res = 0;
    end procedure trr_read_frame

    module procedure trr_skip_frame
        implicit none
        integer(int64) :: position, offset
        logical :: is_open = .false.
        integer :: ierr
        res = -1
        
        if(.not. this%read_header()) return !read header
        
        inquire(this%file, opened = is_open, pos = position) !get position
        if(.not. is_open) return
        
        offset = sum(transfer(this%header%sizes, [integer(int32) ::], 10))
        read(this%file, pos = position + offset, iostat = ierr) !skips the data block
        if(ierr .ne. 0) return

        this%frame%frame_number = this%frame%frame_number + 1
        res = 0
    end procedure trr_skip_frame

    module procedure trr_rewind_file
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
    end procedure trr_rewind_file

    module procedure trr_close_file
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
    end procedure trr_close_file

    module procedure trr_is_open
        implicit none
        inquire(this%file, opened = res)
    end procedure trr_is_open

end submodule TRR_READER