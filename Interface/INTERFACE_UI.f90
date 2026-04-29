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
    implicit none
    
    contains
    
    subroutine print_help()
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
    
end module