! module to handle program switches
	
module switches

contains
	! isFlag(arg)
	! arg - string
	! returns .true. if string starts with '-' and is longer than 2 characters
	logical function isFlag(flag)
		character(len=*), intent(IN) :: flag

		isFlag = .false.

		if(flag(1:1) == '-') then
			 isFlag = .true.
		end if 

		if(len(flag) < 2) isFlag = .false.

	end function isFlag

	! getFlag(arg)
	! arg - string
	! if arg is a switch returns second character after '-', else returns '*'
	character function getFlag(arg)
		character(len=*) :: arg

		if(isFlag(arg)) then
			getFlag = lowercase(arg(2:2))
		else
			getFlag = '*'
		end if
	end function getFlag

	! lowercase(arg)
	! arg - character
	! returns lowercase arg
	character function lowercase(input)
		character, intent(IN) :: input
		lowercase = input
		if(ichar(input) >= ichar('A') .and. ichar(input) <= ichar('Z')) then
			lowercase = char(ichar(input) + 32)
		end if
	end function lowercase

	subroutine getSwitchString(counter, option, targetString)
		integer(kind=8), intent(inout) :: counter
		character, intent(in) :: option
		character(len=*), intent(inout) :: targetString
		character(len = 256) :: arg

		counter = counter + 1
		if(counter > iargc()) then
			print*, "not enough arguments for -", option
			stop
		end if
		call getarg(counter, arg)
		if(.not. isFlag(trim(arg))) then
			targetString = trim(arg)
		else
			print*, trim(arg), " is not a valid option for -", option
			stop
		end if
	end subroutine  getSwitchString

	subroutine getSwitchReal(counter, option, targetReal)
		integer(kind=8), intent(inout) :: counter
		character, intent(in) :: option
		real(kind=8), intent(inout) :: targetReal
		character(len = 256) :: arg
		integer :: ierr
		counter = counter + 1
		if(counter > iargc()) then
			print*, "not enough arguments for -", option
			stop
		end if
		call getarg(counter, arg)
		if(.not. isFlag(trim(arg))) then
			read(arg, *, iostat = ierr) targetReal
			if(ierr .ne. 0) then
				print*, trim(arg), " is not a valid option for -", option
				stop
			end if
		else
			print*, trim(arg), " is not a valid option for -", option
			stop
		end if
	end subroutine getSwitchReal

	subroutine getSwitchInteger(counter, option, targetInteger)
		integer(kind=8), intent(inout) :: counter
		character, intent(in) :: option
		integer(kind=1), intent(inout) :: targetInteger
		character(len = 256) :: arg
		integer :: ierr
		counter = counter + 1
		if(counter > iargc()) then
			print*, "not enough arguments for -", option
			stop
		end if
		call getarg(counter, arg)
		if(.not. isFlag(trim(arg))) then
			read(arg, *, iostat = ierr) targetInteger
			if(ierr .ne. 0) then
				print*, trim(arg), " is not a valid option for -", option
				stop
			end if
		else
			print*, trim(arg), " is not a valid option for -", option
			stop
		end if
	end subroutine getSwitchInteger
	
end module