module FRAME_READERS
    use, intrinsic :: iso_fortran_env, only: real64, int64, int32, real32
    implicit none
    
    type :: fr_atom
        real(real64), dimension(3) ::                       position = (/0.0, 0.0, 0.0/)
        real(real64), dimension(3) ::                       velocity = (/0.0, 0.0, 0.0/)
        character(len=5) ::                                 name = ''
    end type fr_atom

    type :: fr_information
        integer(int64) ::                                   frame_number = 0
        integer(int64) ::                                   n_atoms = 0
        logical ::                                          has_velocities = .false.
    end type fr_information
    
    type, abstract :: frame_reader
        
        
    contains
        procedure(open_), deferred ::                       open_file
        procedure(read_), deferred ::                       read_frame
        procedure(skip_), deferred ::                       skip_frame
        procedure(rewind_), deferred ::                     rewind_file
        procedure(close_), deferred ::                      close_file
        procedure(is_open_), deferred ::                    is_open 
        
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
    
    real(real64), parameter ::                              nmpps_to_hartree = 21.876912635 !nm/ps -> a_0*E_h/(reduced planck)
    real(real64), parameter ::                              nm_to_angstrom = 10

    type(fr_atom), dimension(:), allocatable, protected ::  fr_atoms
    type(fr_information), protected ::                      fr_info

    integer, private ::                                     fr_file = 0, fr_file1 = 0
    integer(int64), private ::                              prev_n_atoms = 0                                       
    
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
    
        procedure ::                                        open_file => trr_open_file
        procedure ::                                        read_frame => trr_read_frame
        procedure ::                                        skip_frame => trr_skip_frame
        procedure ::                                        rewind_file => trr_rewind_file
        procedure ::                                        close_file => trr_close_file
        procedure ::                                        is_open => trr_is_open
        
    end type trr_frame_reader

    interface

        module function trr_open_file(this, filename, filename1) result(res)
            import :: trr_frame_reader
            class(trr_frame_reader), intent(inout) ::       this
            character(*), intent(in) ::                     filename
            character(*), intent(in), optional ::           filename1

            integer ::                                      res
        end function
    
        module function trr_read_frame(this) result(res)
            import :: trr_frame_reader
            class(trr_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function trr_read_frame
    
        module function trr_skip_frame(this) result(res)
            import :: trr_frame_reader
            class(trr_frame_reader), intent(inout) ::       this

            integer ::                                      res        
        end function trr_skip_frame
    
        module function trr_rewind_file(this) result(res)
            import :: trr_frame_reader
            class(trr_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function trr_rewind_file
    
        module function trr_close_file(this) result(res)
            import :: trr_frame_reader
            class(trr_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function trr_close_file
    
        module function trr_is_open(this) result(res)
            import :: trr_frame_reader
            class(trr_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function trr_is_open
    
    end interface
    

    
    !GRO FRAME READER
    
    type, extends(frame_reader) :: gro_frame_reader
        
    contains
    
        procedure ::                                        open_file => gro_open_file
        procedure ::                                        read_frame => gro_read_frame
        procedure ::                                        skip_frame => gro_skip_frame
        procedure ::                                        rewind_file => gro_rewind_file
        procedure ::                                        close_file => gro_close_file
        procedure ::                                        is_open => gro_is_open
        
    end type gro_frame_reader

    interface

        module function gro_open_file(this, filename, filename1) result(res)
            import :: gro_frame_reader
            class(gro_frame_reader), intent(inout) ::       this
            character(*), intent(in) ::                     filename
            character(*), intent(in), optional ::           filename1

            integer ::                                      res
        end function
    
        module function gro_read_frame(this) result(res)
            import :: gro_frame_reader
            class(gro_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function gro_read_frame
    
        module function gro_skip_frame(this) result(res)
            import :: gro_frame_reader
            class(gro_frame_reader), intent(inout) ::       this

            integer ::                                      res        
        end function gro_skip_frame
    
        module function gro_rewind_file(this) result(res)
            import :: gro_frame_reader
            class(gro_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function gro_rewind_file
    
        module function gro_close_file(this) result(res)
            import :: gro_frame_reader
            class(gro_frame_reader), intent(inout) ::       this
            
            integer ::                                      res
        end function gro_close_file
    
        module function gro_is_open(this) result(res)
            import :: gro_frame_reader
            class(gro_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function gro_is_open
    
    end interface
    

    
    !XYZ FRAME READER
    
    type, extends(frame_reader) :: xyz_frame_reader
        
    contains
    
        procedure ::                                        open_file => xyz_open_file
        procedure ::                                        read_frame => xyz_read_frame
        procedure ::                                        skip_frame => xyz_skip_frame
        procedure ::                                        rewind_file => xyz_rewind_file
        procedure ::                                        close_file => xyz_close_file
        procedure ::                                        is_open => xyz_is_open
        
    end type xyz_frame_reader

    interface

        module function xyz_open_file(this, filename, filename1) result(res)
            import :: xyz_frame_reader
            class(xyz_frame_reader), intent(inout) ::       this
            character(*), intent(in) ::                     filename
            character(*), intent(in), optional ::           filename1

            integer ::                                      res
        end function
    
        module function xyz_read_frame(this) result(res)
            import :: xyz_frame_reader
            class(xyz_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function xyz_read_frame
    
        module function xyz_skip_frame(this) result(res)
            import :: xyz_frame_reader
            class(xyz_frame_reader), intent(inout) ::       this

            integer ::                                      res        
        end function xyz_skip_frame
    
        module function xyz_rewind_file(this) result(res)
            import :: xyz_frame_reader
            class(xyz_frame_reader), intent(inout) ::       this

            integer ::                                      res
        end function xyz_rewind_file
    
        module function xyz_close_file(this) result(res)
            import :: xyz_frame_reader
            class(xyz_frame_reader), intent(inout) ::       this
            integer ::                                      res
        end function xyz_close_file
    
        module function xyz_is_open(this) result(res)
            import :: xyz_frame_reader
            class(xyz_frame_reader), intent(inout) ::       this

            logical ::                                      res
        end function xyz_is_open
    
    end interface
    
end module FRAME_READERS