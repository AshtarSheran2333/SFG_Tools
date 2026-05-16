module PRETTY_PRINT
    
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

end module