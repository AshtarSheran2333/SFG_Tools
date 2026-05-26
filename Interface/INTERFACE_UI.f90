!-------------------------------------------------------------------------------
!                               INTERFACE_UI
!
!   This is a set of subroutines that handle program -> console interfacing
!   of the Interface program
!   
!   e.g. handling program options, printing help, printing recap...
!
!-------------------------------------------------------------------------------
module INTERFACE_UI
    use iso_fortran_env
    use SWITCHES
    use UTILS_ERROR
    use PRETTY_PRINT
    use FRAME_READERS
    implicit none
    
#include "utils_error_macros.h"

    integer(kind=int8), parameter ::    ui_filetype_trr = 0,&
                                        ui_filetype_gro = 1,&
                                        ui_filetype_xyz = 2,&
                                        ui_filetype_none = -1

    logical, protected ::               ui_vmd_out  = .false.

    integer(kind=int8), protected ::    ui_filetype = ui_filetype_none

    character(len=128), protected ::    ui_filename1 = "",&
                                        ui_filename2 = ""

    contains
    
    subroutine print_help()
        implicit none

        !TODO might have to be updated
        write(output_unit,f_line)	heading(flat_pattern, "Help dialog of Interface program")
        write(output_unit,f_line)	""                          
        write(output_unit,f_line)	heading(flat_pattern, "Mandatory options")
        write(output_unit,f_line)	""
        write(output_unit,f_line)	"-I(input) <file1> *<file2>"
        write(output_unit,f_1tab)		"Specifies the input files. The options are:"
        write(output_unit,f_2tab)			"*.xyz (positions), *.xyz (velocities)"
        write(output_unit,f_2tab)			"*.gro (positions+velocities)"
        write(output_unit,f_2tab)			"*.trr, *.gro (one frame - can be only positions)"
        write(output_unit,f_line)	""
        write(output_unit,f_line)	heading(flat_pattern, "Optional useful options")
        write(output_unit,f_line)	""
        write(output_unit,f_line)   "-H(help)"
        write(output_unit,f_1tab)       "Prints this help dialog."
        write(output_unit,f_line)   ""
        write(output_unit,f_line)   heading(flat_pattern, "Optional debugging options")
        write(output_unit,f_line)   ""
        write(output_unit,f_line)   "WARNING:"
        write(output_unit,f_1tab)       "These options will create HUGE ASCII files."
        write(output_unit,f_line)   ""
        write(output_unit,f_line)   "-V(vmdout)"
        write(output_unit,f_1tab)       "ASCII dump of interface points in .xyz format"
        write(output_unit,f_line)   ""
        write(output_unit,f_line)   heading(wavy_pattern, "For more details see the documentation")
        write(output_unit,f_line)   ""

    end subroutine print_help
    
    ! goes through program options and sets the program logic
    subroutine evaluate_program_options(fr)
        implicit none
        class(frame_reader), allocatable, intent(inout) :: fr
        character(len=256) :: op, arg, arg1
        character :: option
        integer :: i

        call print_art()
        write(output_unit,f_line) heading(wavy_pattern,"Evaluation of the program options")
        write(output_unit,f_line) ""

        if(command_argument_count() < 1) error_stop("The program needs to be called with some arguments, see help -h")

        ! is help switch present?
        do i = 1, command_argument_count()
            call get_command_argument(i, op)
            if(get_flag(trim(op)) == 'h') then
                call print_help()
                stop
            end if
        end do

        i = 1
        do
            call get_command_argument(i, op)
            option = get_flag(trim(op))

            select case(option)
                case('i')
                    write(output_unit,f_line) "-I(input):"
                    call get_switch_string(i, op, arg)
                    ! INPUT GRO
                    if(index(arg, '.gro') .ne. 0) then
                        write(output_unit, f_1tab) ".gro multiple frame file: "//trim(arg)
                        write(output_unit,f_line) ""
                        ui_filetype = ui_filetype_gro
                        ui_filename1 = trim(adjustl(arg))
                        allocate(gro_frame_reader :: fr)
                        error_io_check(fr%open_file(ui_filename1), "unable to open "//ui_filename1)
                    ! INPUT XYZ
                    else if(index(arg, '.xyz') .ne. 0) then
                        call get_switch_string(i, op, arg1)
                        if(index(arg1, '.xyz') .ne. 0) then
                            write(output_unit,f_1tab) ".xyz position file: "//trim(arg)
                            write(output_unit,f_1tab) ".xyz velocity file: "//trim(arg1)
                            write(output_unit,f_line) ""
                            ui_filetype = ui_filetype_xyz
                            ui_filename1 = trim(adjustl(arg))
                            ui_filename2 = trim(adjustl(arg1))
                            allocate(xyz_frame_reader :: fr)
                            error_io_check(fr%open_file(ui_filename1, ui_filename2), "unable to open files "//ui_filename1//" "//ui_filename2)
                        else
                            !TODO just positions are sufficient
                            error_stop("program needs positions and velocities .xyz files")
                        end if
                    ! INPUT TRR
                    else if(index(arg, '.trr') .ne. 0) then
                        call get_switch_string(i, op, arg1)
                        if(index(arg1, '.gro') .ne. 0) then
                            write(output_unit,f_1tab) ".trr trajectory file: "//trim(arg)
                            write(output_unit,f_1tab) ".gro atleast single frame file: "//trim(arg1)
                            write(output_unit,f_line) ""
                            ui_filetype = ui_filetype_trr
                            ui_filename1 = trim(adjustl(arg))
                            ui_filename2 = trim(adjustl(arg1))
                            allocate(trr_frame_reader :: fr)
                            error_io_check(fr%open_file(ui_filename1, ui_filename2), "unable to open files "//ui_filename1//" "//ui_filename2)
                        else
                            error_stop("program needs .trr and .gro file")
                        end if
                    else
                        error_stop("Unknown input file format: "//trim(arg))
                    end if
                    i = i + 1
                case('v')
                    write(output_unit,f_line) "-V(vmdout) selected"
                    write(output_unit,f_line) ""
                    ui_vmd_out = .true.
                    i = i + 1
                case default
                    write(output_unit,f_line) trim(op)//" skipped - invalid option"
                    write(output_unit,f_line) ""
                    i = i + 1
            end select

            ! all the switches were evaluated
            if(i > command_argument_count()) exit
        end do

        ! check if all the mandatory options were selected
        if(ui_filetype == ui_filetype_none) then
            error_stop("input file(s) must be specified, see help -h")
        end if

        write(output_unit,f_line) heading(flat_pattern,"Evaluation of the program options - DONE")
        write(output_unit,f_line) ""

    end subroutine evaluate_program_options
end module
