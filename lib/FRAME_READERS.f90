module FRAME_READERS
    use, intrinsic :: iso_fortran_env, only: real64, int64, int32, real32
    implicit none

    real(real64), parameter ::                              nmpps_to_hartree = 21.876912635 !nm/ps -> a_0*E_h/(reduced planck)
    real(real64), parameter ::                              nm_to_angstrom = 10

    type :: current_frame_type
        real(real64), dimension(:,:), allocatable ::        positions
        real(real64), dimension(:,:), allocatable ::        velocities
        character(len=5), dimension(:), allocatable ::      names

        integer(int64) ::                                   frame_number = 0
        integer(int64) ::                                   n_atoms = 0
        logical ::                                          has_velocities = .false.
    end type current_frame_type
    
    type, abstract :: frame_reader
        
        !todo exposing current frame as public, do not touch it from outside
        type(current_frame_type), public ::                   frame
        integer, private ::                                 file = 0, file1 = 0
        integer(int64), private ::                          prev_n_atoms = 0                                       
    
    contains
        procedure(open_), deferred, public ::                       open_file
        procedure(read_), deferred, public ::                       read_frame
        procedure(skip_), deferred, public ::                       skip_frame
        procedure(rewind_), deferred, public ::                     rewind_file
        procedure(close_), deferred, public ::                      close_file
        procedure(is_open_), deferred, public ::                    is_open 
        
    end type frame_reader
    
    abstract interface
    
        integer function open_(this, filename, filename1)
            import :: frame_reader
            class(frame_reader), intent(inout) ::           this
            character(*), intent(IN) ::                     filename
            character(*), intent(IN), optional ::           filename1
        end function open_
        
        integer function read_(this)
            import :: frame_reader
            class(frame_reader), intent(inout) ::           this
        end function read_

        integer function skip_(this)
            import :: frame_reader
            class(frame_reader), intent(inout)::            this
        end function skip_

        integer function rewind_(this)
            import :: frame_reader
            class(frame_reader), intent(inout)::            this
        end function rewind_
        
        integer function close_(this) result(res)
            import :: frame_reader
            class(frame_reader), intent(inout)::            this
        end function close_
        
        logical function is_open_(this) result(res)
            import :: frame_reader
            class(frame_reader), intent(inout)::            this
        end function is_open_
        
    end interface
    
    !TRR FRAME READER
    
    type trr_frame_header_sizes
        integer(int32) :: ir_size
        integer(int32) :: e_size
        integer(int32) :: box_size
        integer(int32) :: vir_size
        integer(int32) :: pres_size
        integer(int32) :: top_size
        integer(int32) :: sym_size
        integer(int32) :: pos_size
        integer(int32) :: vel_size
        integer(int32) :: forces_size
    end type trr_frame_header_sizes
        
    type trr_frame_header
        type(trr_frame_header_sizes) :: sizes
            
        integer(int32) :: n_atoms
        integer(int32) :: step_number
        integer(int32) :: nre
        real(real32) :: sim_time
        integer(int32) :: lambda
    end type trr_frame_header
    
    type, extends(frame_reader) :: trr_frame_reader
        
        type(trr_frame_header), private :: header
        
    contains
    
        procedure, public :: open_file => trr_open_file
        procedure, public :: read_frame => trr_read_frame
        procedure, public :: skip_frame => trr_skip_frame
        procedure, public :: rewind_file => trr_rewind_file
        procedure, public :: close_file => trr_close_file
        procedure, public :: is_open => trr_is_open

        procedure, private :: read_header => trr_read_header
        procedure, private :: read_ir => trr_read_ir
        procedure, private :: read_e => trr_read_e
        procedure, private :: read_box => trr_read_box
        procedure, private :: read_vir => trr_read_vir
        procedure, private :: read_pres => trr_read_pres
        procedure, private :: read_top => trr_read_top
        procedure, private :: read_sym => trr_read_sym
        procedure, private :: read_positions => trr_read_positions
        procedure, private :: read_velocities => trr_read_velocities
        procedure, private :: read_forces => trr_read_forces
        procedure, private :: fill_atom_names => trr_fill_atom_names
        
    end type trr_frame_reader

    interface

        module function trr_open_file(this, filename, filename1) result(res)
            class(trr_frame_reader), intent(inout) ::       this
            character(*), intent(in) ::                     filename
            character(*), intent(in), optional ::           filename1

            integer ::                                      res
        end function
    
        module function trr_read_frame(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function trr_read_frame
    
        module function trr_skip_frame(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            integer ::                                      res        
        end function trr_skip_frame
    
        module function trr_rewind_file(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function trr_rewind_file
    
        module function trr_close_file(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function trr_close_file
    
        module function trr_is_open(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_is_open

        module function trr_read_header(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_read_header

        module function trr_read_ir(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_read_ir

        module function trr_read_e(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_read_e

        module function trr_read_box(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_read_box

        module function trr_read_vir(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_read_vir

        module function trr_read_pres(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_read_pres

        module function trr_read_top(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_read_top

        module function trr_read_sym(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_read_sym

        module function trr_read_positions(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_read_positions

        module function trr_read_velocities(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_read_velocities

        module function trr_read_forces(this) result(res)
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_read_forces

        module function trr_fill_atom_names(this, filename) result(res)
            class(trr_frame_reader), intent(inout) ::       this
            character(*) :: filename

            logical ::                                      res
        end function trr_fill_atom_names

    end interface
    

    
    !GRO FRAME READER
    
    type, extends(frame_reader) :: gro_frame_reader
        
    contains
    
        procedure, public :: open_file => gro_open_file
        procedure, public :: read_frame => gro_read_frame
        procedure, public :: skip_frame => gro_skip_frame
        procedure, public :: rewind_file => gro_rewind_file
        procedure, public :: close_file => gro_close_file
        procedure, public :: is_open => gro_is_open
        
        procedure, private :: read_header => gro_read_header

    end type gro_frame_reader

    interface

        module function gro_open_file(this, filename, filename1) result(res)
            class(gro_frame_reader), intent(inout) ::       this
            character(*), intent(in) ::                     filename
            character(*), intent(in), optional ::           filename1

            integer ::                                      res
        end function
    
        module function gro_read_frame(this) result(res)
            class(gro_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function gro_read_frame
    
        module function gro_skip_frame(this) result(res)
            class(gro_frame_reader), intent(inout) ::       this

            integer ::                                      res        
        end function gro_skip_frame
    
        module function gro_rewind_file(this) result(res)
            class(gro_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function gro_rewind_file
    
        module function gro_close_file(this) result(res)
            class(gro_frame_reader), intent(inout) ::       this
            
            integer ::                                      res
        end function gro_close_file
    
        module function gro_is_open(this) result(res)
            class(gro_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function gro_is_open

        module function gro_read_header(this) result(res)
            class(gro_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function gro_read_header

    end interface
    

    
    !XYZ FRAME READER
    
    type, extends(frame_reader) :: xyz_frame_reader
        
    contains
    
        procedure, public :: open_file => xyz_open_file
        procedure, public :: read_frame => xyz_read_frame
        procedure, public :: skip_frame => xyz_skip_frame
        procedure, public :: rewind_file => xyz_rewind_file
        procedure, public :: close_file => xyz_close_file
        procedure, public :: is_open => xyz_is_open
        
        procedure, private :: read_header_with_velocities => xyz_read_header_with_velocities 
        procedure, private :: read_header_no_velocities => xyz_read_header_no_velocities
        procedure, private :: open_with_velocities => xyz_open_with_velocities
        procedure, private :: open_no_velocities => xyz_open_no_velocities

    end type xyz_frame_reader

    interface

        module function xyz_open_file(this, filename, filename1) result(res)
            class(xyz_frame_reader), intent(inout) ::       this
            character(*), intent(in) ::                     filename
            character(*), intent(in), optional ::           filename1

            integer ::                                      res
        end function
    
        module function xyz_read_frame(this) result(res)
            class(xyz_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function xyz_read_frame
    
        module function xyz_skip_frame(this) result(res)
            class(xyz_frame_reader), intent(inout) ::       this

            integer ::                                      res        
        end function xyz_skip_frame
    
        module function xyz_rewind_file(this) result(res)
            class(xyz_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function xyz_rewind_file
    
        module function xyz_close_file(this) result(res)
            class(xyz_frame_reader), intent(inout) ::       this
            integer ::                                      res
        end function xyz_close_file
    
        module function xyz_is_open(this) result(res)
            class(xyz_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function xyz_is_open
    
        module function xyz_read_header_with_velocities(this) result(res)
            class(xyz_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function xyz_read_header_with_velocities

        module function xyz_read_header_no_velocities(this) result(res)
            class(xyz_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function xyz_read_header_no_velocities

        module function xyz_open_with_velocities(this, posfile, velfile) result(res)
            class(xyz_frame_reader), intent(inout) ::       this
            character(*), intent(in) :: posfile, velfile

            integer ::                                      res
        end function xyz_open_with_velocities

        module function xyz_open_no_velocities(this, posfile) result(res)
            class(xyz_frame_reader), intent(inout) ::       this
            character(*), intent(in) :: posfile

            integer ::                                      res
        end function xyz_open_no_velocities

    end interface
    
end module FRAME_READERS
