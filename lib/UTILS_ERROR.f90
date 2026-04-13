module UTILS_ERROR
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    contains

    !stops the program from execution, and prints the message
    subroutine utils_error_io_check(ierr, msg, file, line)
        integer, intent(in) :: ierr
        character(*), intent(in) :: msg
        character(*), intent(in) :: file
        integer, intent(in) :: line
        
        if(ierr /= 0) then
            write(error_unit, "( 'IO ERROR (', I0, '), file: ', A, ' (', I0, '): ', A )") ierr, file, line, msg
            stop
        end if
    end subroutine utils_error_io_check
    
    subroutine utils_error_stop(msg, file, line)
        character(*), intent(in) :: msg
        character(*), intent(in) :: file
        integer, intent(in) :: line

        
        write(error_unit, "( 'ERROR, file: ', A, ' (', I0, '): ', A )") file, line, msg
        stop
    end subroutine utils_error_stop
    
end module UTILS_ERROR
