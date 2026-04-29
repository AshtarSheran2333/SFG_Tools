! module to handle program switches
	
module SWITCHES
	use UTILS_ERROR
	implicit none
    
#include "utils_error_macros.h"

contains
	!TODO private funcion
	! is_flag(arg)
	! arg - string
	! returns .true. if string starts with '-' and is longer than 2 characters
	logical function is_flag(flag)
        implicit none
		character(len=*), intent(IN) :: flag

		is_flag = .false.

		if(flag(1:1) == '-') is_flag = .true.
		if(len(flag) < 2) is_flag = .false.
	end function is_flag

	! get_flag(arg)
	! arg - string
	! if arg is a switch returns second character after '-', else returns '*'
	character function get_flag(arg)
        implicit none
		character(len=*) :: arg

		if(is_flag(arg)) then
			get_flag = lowercase(arg(2:2))
		else
			get_flag = '*'
		end if
	end function get_flag

	! lowercase(arg)
	! arg - character
	! returns lowercase arg
	character function lowercase(input)
        implicit none
		character, intent(IN) :: input
		lowercase = input
		if(ichar(input) >= ichar('A') .and. ichar(input) <= ichar('Z')) then
			lowercase = char(ichar(input) + 32)
		end if
	end function lowercase

	subroutine get_switch_string(counter, option, targetString)
        implicit none
		integer(kind=8), intent(inout) :: counter
		character, intent(in) :: option
		character(len=*), intent(inout) :: targetString
		character(len = 256) :: arg

		counter = counter + 1
		if(counter > iargc()) error_stop("not enough arguments for -"//option)
		call getarg(counter, arg)
		if(.not. is_flag(trim(arg))) then
			targetString = trim(arg)
		else
			error_stop(" is not a valid option for -"//option)
		end if
	end subroutine  get_switch_string

	subroutine get_switch_real64(counter, option, targetReal)
        implicit none
		integer(kind=8), intent(inout) :: counter
		character, intent(in) :: option
		real(kind=8), intent(inout) :: targetReal
		character(len = 256) :: arg
		integer :: ierr
		counter = counter + 1
		if(counter > iargc()) error_stop("not enough arguments for -"//option)
		call getarg(counter, arg)
		if(.not. is_flag(trim(arg))) then
			read(arg, *, iostat = ierr) targetReal
			if(ierr .ne. 0) then
				error_stop(" is not a valid option for -"//option)
			end if
		else
			print*, trim(arg), " is not a valid option for -", option
			stop
		end if
	end subroutine get_switch_real64

	subroutine get_switch_int8(counter, option, targetInteger)
        implicit none
		integer(kind=8), intent(inout) :: counter
		character, intent(in) :: option
		integer(kind=1), intent(inout) :: targetInteger
		character(len = 256) :: arg
		integer :: ierr
		counter = counter + 1
		if(counter > iargc()) then
			error_stop("not enough arguments for -"//option)
		end if
		call getarg(counter, arg)
		if(.not. is_flag(trim(arg))) then
			read(arg, *, iostat = ierr) targetInteger
			if(ierr .ne. 0) then
				print*, trim(arg), " is not a valid option for -", option
				stop
			end if
		else
			print*, trim(arg), " is not a valid option for -", option
			stop
		end if
	end subroutine get_switch_int8

	subroutine get_switch_int64(counter, option, targetInteger)
        implicit none
		integer(kind=8), intent(inout) :: counter
		character, intent(in) :: option
		integer(kind=8), intent(inout) :: targetInteger
		character(len = 256) :: arg
		integer :: ierr
		counter = counter + 1
		if(counter > iargc()) then
			error_stop("not enough arguments for -"//option)
			stop
		end if
		call getarg(counter, arg)
		if(.not. is_flag(trim(arg))) then
			read(arg, *, iostat = ierr) targetInteger
			if(ierr .ne. 0) then
				print*, trim(arg), " is not a valid option for -", option
				stop
			end if
		else
			print*, trim(arg), " is not a valid option for -", option
			stop
		end if
	end subroutine get_switch_int64
	
end module

