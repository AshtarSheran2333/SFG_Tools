!this is extremely ugly, stitched fast, it works... can be refactored later
program Corr_to_spectrum

	use bdata
	
	implicit none
	
	type(boxdata) ::                                            bd
	
	real*8, allocatable, dimension(:) ::                        corr
	
	real*8 ::                                                   arbitrary_real,&
																filterpower
	
	integer*8 ::                                                i,&
																j,&
																k,&
																corr_file_len,&
																maxlag=0
	
	integer ::                                                  ierr,&
																input_file_unit,&
																output_file_unit
	
	character(len=256) ::                                       line
	character(len=256) ::										output=""
	
	real*8, parameter ::                                        pi = dacos(-1.d0),&
																kb = 1.386658d-23,&
																c = 2.99792458d10
	
	

!________________________BODY OF THE PROGRAM______________________
	
	call evaluate_program_switches()
	
	call bd%read_boxdata()
	
	allocate(corr(maxlag))
	
	do while(ierr .ne. 0)
		read(input_file_unit,"(A)",iostat = ierr) line
		if(index(adjustl(line),"#") == 1) then
			cycle
		else
			read(line,*), j, corr(i)
		end if
	end do

	i = 1
	do while(ierr .eq. 0)
			read(input_file_unit,"(A)",iostat = ierr) line
			if(ierr .ne. 0) then
				exit
			else
				if(index(adjustl(line),"#") == 1) then
					cycle
				else
					read(line,*,iostat = ierr), arbitrary_real, corr(i)
					if(ierr .ne. 0) then
						exit
					else
						i = i + 1
					end if
				end if
			end if
		end do

	if(i-1 == maxlag) then
		print*, "Reading of the input file succesful"
	else
		print*, "ERROR: inconsitency in data reading" 
		stop
	end if
		
	if(bd%filter <= 0) then
		!filter <= 0 makes no sense -> Simones default version...
		!default filter value by Simone after the change of filter
		filterpower = maxlag/(10*sqrt(2.0))
	else
		!filter is set according to the value in boxdata (can be overriden by -filter switch in program switches evaluation) 
		filterpower = 1000*bd%filter/bd%DT
	end if
	
	call make_spectrum(corr, trim(output))
	
	print*, "done"
contains
	
subroutine make_spectrum(corr, name)
	real(8), dimension(:), intent(IN) ::                    corr
	
	character(len=*), intent(IN) ::                         name
	
	complex(8), allocatable, dimension(:) ::                spectrum
	
	complex(8) ::                                           eiomegat,&
															eiomegat1
	
	real(8) ::                                              time,&
															wavenumber,&
															omega,&
															beta
	
	integer(8) ::                                           f,&
															f_spectrum,&
															f_shiftedspectrum, ierr
	
	f = bd%FREQ/bd%DFREQ

	allocate(spectrum(0:f))
	spectrum = 0
	beta = 1/(kb*bd%TEMPERATURE)
	
	do i = 0, f
		wavenumber = dble(i)*bd%DFREQ !cm^-1
		omega = 2.d0*pi*c*wavenumber !s^-1
		
		!LEFT RULE
		!do j = 1, maxlag
		!    time=bd%DT*dble(j-1)*1.d-15 !s
		!    eiomegat = dcmplx(dcos(dble(omega*time)), dsin(dble(omega*time)))
		!    spectrum(i) = spectrum(i) + eiomegat * corr(j) * filter(j, filterpower) * bd%DT
		!end do
		
		!TRAPEZOIDAL RULE - why not?!?!
		do j = 1, maxlag
			time=bd%DT*dble(j-1)*1.d-15 !s
			eiomegat = dcmplx(dcos(dble(omega*time)), dsin(dble(omega*time)))
			time=bd%DT*dble(j)*1.d-15 !s
			eiomegat1 = dcmplx(dcos(dble(omega*time)), dsin(dble(omega*time)))
			if(j < maxlag) then
				spectrum(i) = spectrum(i) + 0.5d0 * ( eiomegat * corr(j) * filter(j, filterpower) + eiomegat1 * corr(j+1) * filter(j+1, filterpower) ) * bd%DT
			else !add one artificial point to be 0 for the trap rule...
				spectrum(i) = spectrum(i) + 0.5d0 * eiomegat * corr(j) * filter(j, filterpower) * bd%DT
			end if
		end do
		
		spectrum(i) = dcmplx(0,-1) * beta * spectrum(i) / omega
	end do
	
	print*,""
	print*, trim(name)//"-spectrum.dat"
	open(newunit = f_spectrum, file=trim(name)//"-spectrum.dat", recl=120, iostat = ierr)
	if(ierr .ne. 0) then
		print"(a,a,a)", "ERROR: unable to create file '", trim(name)//"-spectrum.dat", "'"
		stop
	end if
	write(f_spectrum,"(A,F8.3)"), "#FILTER PARAMETER: ", (filterpower/1000)*bd%DT !in fs
	write(f_spectrum,"(A)") "#wavelength cm-1, re, im, abs, phase deg"
	print*, trim(name)//"-shiftedspectrum.dat"
	print*, trim(name)//"-shiftedspectrum.dat is shifted by:"
	print*, "(Re, Im)", -spectrum(f)
	print*, "which is value of the spectrum at frequency ", f, "cm^-1"
	open(newunit = f_shiftedspectrum, file=trim(name)//"-shiftedspectrum.dat", recl=120, iostat = ierr)
	if(ierr .ne. 0) then
		print"(a,a,a)", "ERROR: unable to create file '", trim(name)//"-shiftedspectrum.dat", "'"
		stop
	end if
	write(f_shiftedspectrum,"(A,F8.3)"), "#FILTER PARAMETER: ", (filterpower/1000)*bd%DT !in fs
	write(f_shiftedspectrum,"(A)") "#wavelength cm-1, re, im, abs, phase deg"
	
	do i=0,f
		!according to Khatib equation (3) gets the second order susceptibility
		wavenumber=dble(i)*bd%DFREQ       !cm-1
	
		write(f_spectrum,*) int(wavenumber),&
					real(spectrum(i)),&
					imag(spectrum(i)),&
					abs(spectrum(i)),&
					datan(imag(spectrum(i))/real(spectrum(i)))*180/pi
	 
		write(f_shiftedspectrum,*) wavenumber,&
					real(spectrum(i)-spectrum(f)),&
					imag(spectrum(i)-spectrum(f)),&
					abs(spectrum(i)-spectrum(f)),&
					datan(imag(spectrum(i)-spectrum(f))/real(spectrum(i)-spectrum(f)))*180/pi
	enddo
	
	close(f_spectrum)
	close(f_shiftedspectrum)
	
	deallocate(spectrum)
end subroutine make_spectrum    
	
real*8 function filter(time,tau)

	integer*8, intent(IN) ::                                time
	real*8, intent(IN) ::                                   tau
	
	filter = dexp(-(dble(time)/dble(tau))**2)
end function filter

subroutine evaluate_program_switches()

	character(len=60), allocatable ::                           arg,&
																arg1
	allocate(arg)
	allocate(arg1)

	i = -1
	do
		!if there is a help just skip to the help dialog
		if(i < 0) then
			do j=1, iargc()
				if(arg .ne. "-h") then
					call getarg(j, arg)
				end if
			end do
			i = 0
			if(arg .ne. "-h") arg = ""
		else !else loop through all the options
			i = i+1
			call getarg(i, arg)
		end if
		
	select case(arg)
	case("-i")
		i = i+1
		call getarg(i, arg)
		
		open(newunit = input_file_unit, file = arg, status = 'old', iostat = ierr)
		if(ierr .ne. 0) then
			print "(a,a,a)", "Input file ", trim(arg), " does not exist"
			stop
		end if
		
		!get number of lines and skip 1st line
		
		do while(ierr .eq. 0)
			read(input_file_unit,"(A)",iostat = ierr) line
			if(ierr .ne. 0) then
				exit
			else
				if(index(adjustl(line),"#") == 1) then
					cycle
				else
					read(line,*,iostat = ierr), arbitrary_real, arbitrary_real
					if(ierr .ne. 0) then
						exit
					else
						maxlag = maxlag + 1
					end if
				end if
			end if
		end do
		
		if(maxlag == 0) then
			print*, "Not enough datapoints in the input file"
			stop
		end if
		rewind(input_file_unit)
		read(input_file_unit,*)
		
	case("-o")
		i = i+1
		call getarg(i, arg)
		output = arg
		
	case("-filter")
		i = i+1
		call getarg(i, arg)
		read(arg,*, iostat = ierr) bd%filter
		if(ierr .ne. 0) then
			print*, "-filter <arg> does not have a valid argument"
		end if
	
	case("-h")
		
	case("")
		
	case default 
		print*, "program needs atleast -i -*filename*"
		stop
	end select
	
	if(i >= iargc()) then
		exit
	end if
end do

deallocate(arg)
deallocate(arg1)
print*, "_______________________________________________________________________________"    

end subroutine evaluate_program_switches

end program Corr_to_spectrum

