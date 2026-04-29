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
    implicit none
    
#include "utils_error_macros.h"

    integer(kind=int8), parameter ::    ui_filetype_trr = 0,&
                                        ui_filetype_gro = 1,&
                                        ui_filetype_xyz = 2,&
                                        ui_filetype_none = -1

    logical, protected ::   ui_no_check = .false.,&
                            ui_vmd_out  = .false.

    integer(kind=int8), protected ::    ui_filetype = ui_filetype_none

    character(len=128), protected ::    ui_filename1 = "",&
                                        ui_filename2 = ""

    contains
    
    subroutine print_help()
        implicit none
        character(len=*), parameter ::  f_line = "(A)",&
                                        f_1tab = "(T4, A)",&
                                        f_2tab = "(T8,A)"

        !TODO might have to be updated
        write(output_unit,f_line)	""
        write(output_unit,f_line)	"__________________________Help dialog of Interface_____________________________"
        write(output_unit,f_line)	""                          
        write(output_unit,f_line)	"_____________________________Mandatory options_________________________________"
        write(output_unit,f_line)	""
        write(output_unit,f_line)	"-I(input) <file1> *<file2>"
        write(output_unit,f_1tab)		"specifies the input files. Options are:"
        write(output_unit,f_2tab)			"*.xyz (positions), *.xyz (velocities)"
        write(output_unit,f_2tab)			"*.gro (positions+velocities)"
        write(output_unit,f_2tab)			"*.trr, *.gro (one frame - can be only positions)"
        write(output_unit,f_line)	""
        write(output_unit,f_line)	"___________________________Optional useful options_____________________________"
        write(output_unit,f_line)	""
        write(output_unit,f_line)   "-G(good)"
        write(output_unit,f_1tab)       "skips checking of the input files"
        write(output_unit,f_line)   ""
        write(output_unit,f_line)   "-H(help)"
        write(output_unit,f_1tab)       "prints this help dialog"
        write(output_unit,f_line)   ""
        write(output_unit,f_line)   "__________________________Optional debugging options___________________________"
        write(output_unit,f_line)   ""
        write(output_unit,f_line)   "WARNING:"
        write(output_unit,f_1tab)       "These options will create HUGE ASCII files."
        write(output_unit,f_1tab)       "These options are kept in the code mainly for historical reasons/debugging."
        write(output_unit,f_line)   ""
        write(output_unit,f_line)   "-V(vmdout)"
        write(output_unit,f_1tab)       "ASCII dump of interface points in .xyz format"
        write(output_unit,f_line)   ""
        write(output_unit,f_line)   "for more details see documentation"
    end subroutine print_help
    
    ! goes through program options and sets the program logic
    ! frame reader as an argument...
    subroutine evaluate_program_options()
        implicit none
        character(len=256) :: op, arg, arg1
        character :: option
        integer(kind=int64) :: i

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
                    call get_switch_string(i, op, arg)
                    ! INPUT GRO
                    if(index(arg, '.gro') .ne. 0) then
                        print*, "input file:"
                        print*, ""
                        print"(tr5,a,tr5,a)", ".gro multiple frame file: ", trim(arg)
                        print*, ""
                        ui_filetype = ui_filetype_gro
                        ui_filename1 = trim(adjustl(arg))
                    ! INPUT XYZ
                    else if(index(arg, '.xyz') .ne. 0) then
                        call get_switch_string(i, op, arg1)
                        if(index(arg1, '.xyz') .ne. 0) then
                            print*, "input files:"
                            print*, ""
                            print"(tr5,a,tr5,a)", ".xyz position file: ", trim(arg)
                            print"(tr5,a,tr5,a)", ".xyz velocity file: ", trim(arg1)
                            print*, ""
                            ui_filetype = ui_filetype_xyz
                            ui_filename1 = trim(adjustl(arg))
                            ui_filename2 = trim(adjustl(arg1))
                        else
                            !TODO just positions are sufficient
                            error_stop("program needs position and velocity .xyz files in order position file velocity file")
                        end if
                    ! INPUT TRR
                    else if(index(arg, '.trr') .ne. 0) then
                        call get_switch_string(i, op, arg1)
                        if(index(arg1, '.gro') .ne. 0) then
                            print*, "input files:"
                            print*, ""
                            print"(tr5,a,tr5,a)", ".trr trajectory file: ", trim(arg)
                            print"(tr5,a,tr5,a)", ".gro atleast single frame file: ", trim(arg1)
                            print*, ""
                            ui_filetype = ui_filetype_trr
                            ui_filename1 = trim(adjustl(arg))
                            ui_filename2 = trim(adjustl(arg1))
                        else
                            error_stop("program needs .trr and .gro files in order .trr file .gro file")
                        end if
                    else
                        error_stop("Unknown input file format: "//trim(arg))
                    end if

                    i = i + 1
                case('g')
                    ui_no_check = .true.
                    i = i + 1
                case('v')
                    ui_no_check = .true.
                    i = i + 1
                case default
                    print*, "skipping invalid option ", trim(op)
                    i = i + 1
            end select

            ! all the switches were evaluated
            if(i > command_argument_count()) exit
        end do

        print*, ""

        ! check if all the mandatory options were selected
        if(ui_filetype == ui_filetype_none) then
            error_stop("input file(s) must be specified, see help -h")
        end if

    end subroutine evaluate_program_options
end module