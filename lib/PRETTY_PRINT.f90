module PRETTY_PRINT
    use iso_fortran_env
    implicit none
    
    character(len=*), parameter ::  f_line = "(A)",&
                                    f_1tab = "(T4, A)",&
                                    f_2tab = "(T8,A)"

    character(len=80), parameter :: &
        wavy_pattern = "~._.~'^'~._.~'^'~._.~'^'~._.~'^'~._.~'^'~._.~'^'~._.~'^'~._.~'^'~._.~'^'~._.~'^'",&
        flat_pattern = "________________________________________________________________________________"

    contains
    
    subroutine print_art()
        write(output_unit,f_line) ":----=+++==-::...                      #### #####  ###                          "
        write(output_unit,f_line) ".::---=++++=--::..                    #     #     #                             "
        write(output_unit,f_line) "..:::---=++++=-::...                   ###  ####  #  ##                         " 
        write(output_unit,f_line) "...:::--==+++=--::...                     # #     #   #                         "
        write(output_unit,f_line) " ....:::--=++++=--::..                ####  #      ###                       ..."
        write(output_unit,f_line) "   ....::---=+++==-::...                                                   ....."
        write(output_unit,f_line) "     ....::--==+++=--::...     *****  ***   ***  *      ****            ......:-"
        write(output_unit,f_line) "      ....:::--=++++=-:::..      *   *   * *   * *     *             ......:::-="
        write(output_unit,f_line) "....    ....:::--=+++==-::...    *   *   * *   * *      ***       ......:::---=+"
        write(output_unit,f_line) "-:...... .....:::-=++++=--::..   *   *   * *   * *         *    .....:::--==+++="
        write(output_unit,f_line) "--::::..........::--=+++==-:::.  *    ***   ***  ***** **** ......:::--==++++==:"
        write(output_unit,f_line) "+----:::::........::--=+++=--::...                        .....:::--==++++==---:"
        write(output_unit,f_line) "+++===---::::.......::-==++==--::..                    .....:::--==++++==---:::."
        write(output_unit,f_line) "-==+++++==---::::....:::-=+++=--::...                ....:::--=+++++==---:::... "
        write(output_unit,f_line) ":::--==+++++===--::::..:::-=+++=--::...           ....:::--=++++==---::::....   "
        write(output_unit,f_line) " ...:::--==++++++==--::::::-==++==--::::...........:::--=++++==---:::....       "
        write(output_unit,f_line) "    ....:::--==++++++==--::::-=+++===----::::::::::-==++++==--::::....          "
        write(output_unit,f_line) "        .....:::--=++++++==--:--=+++======-------==++++=---:::.....             "
        write(output_unit,f_line) "             ....:::--==++++++====+++++++======++++==--::::....                 "
        write(output_unit,f_line) "                  ....::--==++++++++++++++++++++==--::::....                    "
        write(output_unit,f_line) "                      ....::-==+++++++++++++++==--:::....                       "
        write(output_unit,f_line) "                         ..::--==++++++++++++===--::..                          "
    end subroutine print_art

    
    function heading(pattern, str) result(res)
        character(*), intent(in) :: str
        character(len=80), intent(in) :: pattern
        character(len=80) :: res
        integer :: pos, l

        res = pattern

        l = min(len(str), 68)
        if(l == 0) return
        pos = 40 - l/2 - mod(l,2) + 1
        res(pos-1:pos-1) = " "
        res(pos:pos+l-1) = str(1:l)
        res(pos+l:pos+l) = " "
    end function heading

    subroutine print_main_loop_progress(iteration, steps)
        integer(int64), intent(in) :: iteration, steps
        integer(int64) :: current_value
        real(real32), save :: last_value = -1

        current_value = int(1000 * real(iteration) / real(steps))
        
        if(current_value .ne. last_value) then
            write(output_unit,"(A,' ',F5.1,'%')") "Progress:", current_value / 10.0
            last_value = current_value
        end if
    end subroutine

end module