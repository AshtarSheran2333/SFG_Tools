submodule (FRAME_READERS) TRR_READER
    use, intrinsic :: iso_fortran_env, only: output_unit
    use UTILS_ERROR
    implicit none

#include "utils_error_macros.h"
    
    contains
    
    function read_header(reader) result(res)
        type(trr_frame_reader), intent(inout) :: reader
        logical :: res
        integer(int32) :: skip = 0
        character(:), allocatable :: cw
        logical :: is_open
        integer :: ierr
        
        res = .false.

        inquire(fr_file, opened = is_open)

        if(.not. is_open) return
        
        !magic number
        read(fr_file, iostat = ierr) skip !1993
        if(ierr .ne. 0) return

        if(skip /= 1993) then
            write(output_unit,*) "corrupted trr file header"
            return
        end if

        read(fr_file, iostat = ierr) skip !max str len
        if(ierr .ne. 0) return

        read(fr_file, iostat = ierr) skip !actual str len
        if(ierr .ne. 0) return

        allocate(character(len=skip) :: cw)
        
        read(fr_file, iostat = ierr) cw !"GMX_trn_file"
        if(ierr .ne. 0) return

        if(cw .ne. "GMX_trn_file") then
            write(output_unit,*) "unexpected trr string"
            return
        end if

        deallocate(cw)
        
        !àead the whole header in one go
        associate(h => reader%header)
        read(fr_file, iostat = ierr) &
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
    end function
    
    function read_ir(reader) result(res) !not implemented - just skipping the data
        class(trr_frame_reader) :: reader
        logical :: res
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(reader%header%sizes%ir_size .eq. 0) then 
            res = .true.
            return !nothing to read
        end if

        inquire(fr_file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(fr_file, pos = position + reader%header%sizes%ir_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong
        
        res = .true.
    end function read_ir

    function read_e(reader) result(res) !not implemented - just skipping the data
        class(trr_frame_reader) :: reader
        logical :: res
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(reader%header%sizes%e_size .eq. 0) then 
            res = .true.
            return !nothing to read
        end if
        
        inquire(fr_file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(fr_file, pos = position + reader%header%sizes%e_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong
        
        res = .true.
    end function read_e

    function read_box(reader) result(res) !not implemented - just skipping the data
        class(trr_frame_reader) :: reader
        logical :: res
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(reader%header%sizes%box_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(fr_file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(fr_file, pos = position + reader%header%sizes%box_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end function read_box

    function read_vir(reader) result(res) !not implemented - just skipping the data
        class(trr_frame_reader) :: reader
        logical :: res
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(reader%header%sizes%vir_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(fr_file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(fr_file, pos = position + reader%header%sizes%vir_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end function read_vir

    function read_pres(reader) result(res) !not implemented - just skipping the data
        class(trr_frame_reader) :: reader
        logical :: res
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(reader%header%sizes%pres_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(fr_file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(fr_file, pos = position + reader%header%sizes%pres_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end function read_pres

    function read_top(reader) result(res) !not implemented - just skipping the data
        class(trr_frame_reader) :: reader
        logical :: res
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(reader%header%sizes%top_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(fr_file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(fr_file, pos = position + reader%header%sizes%top_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end function read_top
    
    function read_sym(reader) result(res) !not implemented - just skipping the data
        class(trr_frame_reader) :: reader
        logical :: res
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(reader%header%sizes%sym_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(fr_file, pos = position, opened = is_open)
        if(.not. is_open) return
        
        read(fr_file, pos = position + reader%header%sizes%sym_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end function read_sym
    
    function read_positions(reader) result(res)
        class(trr_frame_reader) :: reader
        logical :: res
        integer :: ierr
        real(real32), dimension(3) :: read_sp
        integer(int32) :: element_sz
        integer(int32) :: i
        logical :: is_open
        
        res = .false.
        
        if(reader%header%sizes%pos_size .eq. 0) then
            return !nothing to read
        end if

        inquire(fr_file, opened = is_open)
        if(.not. is_open) return
        
        element_sz = reader%header%sizes%pos_size / (3 * reader%header%n_atoms)
        
        if( (element_sz * 3 * reader%header%n_atoms) .ne. reader%header%sizes%pos_size ) &
            error_stop("position size is not congruent with the number of atoms")
        
        if(reader%header%n_atoms .ne. size(fr_atoms)) &
            error_stop("trr reader internal n_atoms does not match the trr header n_atoms")
        
        if( (element_sz .ne. 4) .and. (element_sz .ne. 8) ) &
            error_stop("trr positions are neither single nor double precision")
        
        if(element_sz .eq. 8) then
            !read double precision
            read(fr_file, iostat = ierr) (fr_atoms(i)%position, i = 1, reader%header%n_atoms) 
        else
            !read single precision, store in double precision
            do i = 1, reader%header%n_atoms
                  read(fr_file, iostat = ierr) read_sp
                  if(ierr .ne. 0) exit
                  fr_atoms(i)%position = read_sp
            end do
        end if

        if(ierr .ne. 0) return !something went wrong

        !convert the positions
        do i = 1, size(fr_atoms)
            fr_atoms(i)%position = fr_atoms(i)%position * nm_to_angstrom
        end do
        
        res = .true.
    end function read_positions
    
    function read_velocities(reader) result(res)
        class(trr_frame_reader) :: reader
        logical :: res
        integer :: ierr
        real(real32), dimension(3) :: read_sp
        integer(int32) :: element_sz
        integer(int32) :: i
        logical :: is_open
        
        res = .false.
        
        if(reader%header%sizes%vel_size .eq. 0) then
            if(fr_info%has_velocities) return !error, config changed
            res = .true.
            fr_info%has_velocities = .false.
            return !nothing to read
        end if

        inquire(fr_file, opened = is_open)
        if(.not. is_open) return
        
        element_sz = reader%header%sizes%vel_size / (3 * reader%header%n_atoms)
        
        if( (element_sz * 3 * reader%header%n_atoms) .ne. reader%header%sizes%vel_size ) &
            error_stop("trr velocity size is not congruent with the number of atoms")
        
        if(reader%header%n_atoms .ne. size(fr_atoms)) &
            error_stop("trr reader internal n_atoms does not match the trr header n_atoms")
        
        if( (element_sz .ne. 4) .and. (element_sz .ne. 8) ) &
            error_stop("trr velocities are neither single nor double precision")
        
        if(element_sz .eq. 8) then
            !read double precision
            read(fr_file, iostat = ierr) (fr_atoms(i)%velocity, i = 1, reader%header%n_atoms) 
        else
            !read single precision, store in double precision
            do i = 1, reader%header%n_atoms
                  read(fr_file, iostat = ierr) read_sp
                  if(ierr .ne. 0) exit
                  fr_atoms(i)%velocity = read_sp
            end do
        end if

        if(ierr .ne. 0) return !something went wrong

        !convert the velocities
        do i = 1, size(fr_atoms)
            fr_atoms(i)%velocity = fr_atoms(i)%velocity * nmpps_to_hartree
        end do

        fr_info%has_velocities = .true.
        res = .true.
    end function read_velocities

    function read_forces(reader) result(res) !not implemented - just skipping the data
        class(trr_frame_reader) :: reader
        logical :: res
        integer :: ierr
        integer(int64) :: position
        logical :: is_open
        
        res = .false.
        
        if(reader%header%sizes%forces_size .eq. 0) then
            res = .true.
            return !nothing to read
        end if 
    
        inquire(fr_file, pos = position, opened = is_open)
        
        if(.not. is_open) return
        
        read(fr_file, pos = position + reader%header%sizes%forces_size, iostat = ierr) !just skip 
        if(ierr .ne. 0) return !something went wrong

        res = .true.
    end function read_forces

    module procedure trr_open_file
        integer :: ierr, file1
        logical :: is_open
        res = -1

        ! not necessary, but I want to have the atomnames!
        if(.not. present(filename1)) error_stop("trr_open_file must be called with both arguments")
        
        if(allocated(fr_atoms)) then
            deallocate(fr_atoms)
        end if
        
        fr_info%frame_number = 0
        fr_info%n_atoms = 0
        fr_info%has_velocities = .false.
        
        inquire(fr_file, opened = is_open)
        if(is_open) then
            close(fr_file)
        end if
        
        write(output_unit,'( "Opening ", A, " file...")') trim(filename)

        open(newunit = fr_file, file = filename, form = "unformatted",&
                access = 'stream', convert = 'big_endian', status = 'old', iostat = ierr)
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        
        if(.not. read_header(this)) return
        
        fr_info%n_atoms = this%header%n_atoms
        prev_n_atoms = fr_info%n_atoms

        if(this%header%sizes%vel_size > 0) fr_info%has_velocities = .true.

        allocate(fr_atoms(this%header%n_atoms), stat = ierr) !allocate space for atoms
        if(ierr .ne. 0) return
        
        rewind(fr_file, iostat = ierr) !rewind
        if(ierr .ne. 0) return

        res = 0
    end procedure trr_open_file

    module procedure trr_read_frame
        logical :: is_open
        
        res = -1
        inquire(fr_file, opened = is_open)
        if(.not. is_open) return

        res = -2; if(.not. read_header(this)) return
        
        !check wheather the number of atoms changed, if yes, reallocate fr_atoms
        if(this%header%n_atoms .ne. prev_n_atoms) then
            error_stop("trr reader does not support variable number of atoms") 
            !well, it does if the line above is commented, but I dont want this functionality
            !since updating the atomnames would be pain in the ass
            if(allocated(fr_atoms)) then
                deallocate(fr_atoms)
            end if
            allocate(fr_atoms(this%header%n_atoms))
            prev_n_atoms = this%header%n_atoms
        end if
        
        res = -3; if(.not. read_ir(this)) return
        res = -4; if(.not. read_e(this)) return
        res = -5; if(.not. read_box(this)) return
        res = -6; if(.not. read_vir(this)) return
        res = -7; if(.not. read_pres(this)) return
        res = -8; if(.not. read_top(this)) return
        res = -9; if(.not. read_sym(this)) return
        res = -10; if(.not. read_positions(this)) return
        res = -11; if(.not. read_velocities(this)) return
        res = -12; if(.not. read_forces(this)) return

        fr_info%frame_number = fr_info%frame_number + 1
        res = 0;
    end procedure trr_read_frame

    module procedure trr_skip_frame
        integer(int64) :: position, offset
        logical :: is_open = .false.
        integer :: ierr
        res = -1
        
        if(.not. read_header(this)) return !read header
        
        inquire(fr_file, opened = is_open, pos = position) !get position
        if(.not. is_open) return
        
        offset = sum(transfer(this%header%sizes, [integer(int32) ::], 10))
        read(fr_file, pos = position + offset, iostat = ierr) !skips the data block
        if(ierr .ne. 0) return

        fr_info%frame_number = fr_info%frame_number + 1
        res = 0
    end procedure trr_skip_frame

    module procedure trr_rewind_file
        logical :: is_open
        integer :: ierr
        
        inquire(fr_file, opened = is_open)
        res = -1; if(.not. is_open) return

        rewind(fr_file, iostat = ierr)
        if(ierr .ne. 0) then
            res = ierr
            return
        end if
        
        fr_info%frame_number = 0
        res = 0
    end procedure trr_rewind_file

    module procedure trr_close_file
        logical :: is_open
        
        inquire(fr_file, opened = is_open)
        res = -1; if(.not. is_open) return
        
        close(fr_file, iostat = res)
        if(res .ne. 0) return

        fr_file = 0

        if(allocated(fr_atoms)) deallocate(fr_atoms)
        
        fr_info%n_atoms = 0
        fr_info%frame_number = 0
        fr_info%has_velocities = .false.
        res = 0 
    end procedure trr_close_file

    module procedure trr_is_open
        inquire(fr_file, opened = res)
    end procedure trr_is_open

end submodule TRR_READER