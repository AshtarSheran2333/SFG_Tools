module bdata

implicit none
	
	type boxdata
		
		!PRIVATE variables
		integer, private ::                                 fileUnit,&
															layersCount = 1

		character(len=3), private ::						polarization = "SSP"	
		
		!PUBLIC variables (boxdata parameters)
		real*8, dimension(3) ::                             box_dimensions = (/-1,-1,-1/),& !Angstrom
															interface_volume_element = (/0.5d0, 0.5d0, 0.25d0/) !Angstrom
		
		real*8, dimension(7) ::								layers_limits = (/0,0,0,0,0,0,0/) !Angstrom
															
		
		real*8 ::                                           FREQ = 4000,& !cm-1
															DFREQ = 1,& !cm-1
															DT = -1,& !fs
															CORRLEN = 0,& !ps
															FILTER = 0,& !ps
															TEMPERATURE = 300,& !K
															!some parameters related to binning and ranges of the plots...
															hbhist_anglestart = -1d0,& !cos angle
															hbhist_angleend = 1d0,& !cos angle
															hbhist_diststart = 2.3,& !min O-O distance A
															hbhist_distend = 3.2,& !max O-O distance A                                                          
															dipole_rstart = -5.d0,& !min distance from instasurf A
															dipole_rend = 80.d0,& !max distance from instasurf A
															dipole_rbinwidth = 0.025d0,& !binwidth of dipole_vs_r
															density_min_r = -4,& !min distance from instasurf A
															density_max_r = 25,& !min distance from instasurf A
															density_tol = 0,& !density bin epsilon A
															interface_density = 0.5 !the amount of water bulk density to be IS -
		
		integer ::                                          NATOM=-1,& !number of atoms -
															NO=-1,& !number of oxygens -
															NSTEP=-1,& !number of steps -
															INTERFACE_SKIP=1,& !number of frames skipped when evaluating IS -
															CROSS_SKIP=1,& !number time beginnings to skip for crossterms -
															SELF_SKIP=1,& !number of time beginnings to skip for selfterms -
															INTERFACE_PUSHBACK=100000000,& !-
															!some parameters related to binning and ranges of the plots...
															hbhist_anglediv = 100,& !number of divisions angleend-anglestart
															hbhist_distdiv = 45,& !number of divisions distend-diststart
															HBHIST_SKIP=1,& !number of frames skipped when evaluating hydrogen bonds -
															density_bin = 5, & !number of bins in 1 A
															!polarization setup variables
															P = 1, & !chi2_{pqr}
															Q = 1, & !chi2_{pqr}
															R = 3 !chi2_{pqr}

		logical ::											cancel_interface_calculation = .false. !interface loads existing binder to bypass the interface calculation to get only the density
															
		character(10) ::                                    hydroxyl_metal = "X"
		
	contains

		procedure ::            read_boxdata
		procedure ::            get_maxlag
		procedure ::            get_layersCount
		procedure, private ::	print_recap
	end type
	
contains

subroutine read_boxdata(this)
	class(boxdata) ::                                       this
	character(len = 255) ::                                 line
	integer ::                                              read_status = 0,&
															ierr,&
															i
	logical ::                                              error = .false.
	
	open(newunit = this%fileUnit, file = "BOXDATA", form = 'formatted',status = 'old', IOSTAT = read_status)
	print*, "Opening BOXDATA..."
	print*, ""
	
	!if opened read one line if not print error and end...
	if(read_status .ne. 0) then
		print*, "Problem Opening BOXDATA file (file probably does not exist)"
		stop
	end if
	
	read(this%fileUnit,'(A)', IOSTAT = read_status) line
	
	!evaluate all the lines in BOXDATA
	do while (read_status .EQ. 0)
			
		select case(line)
		
			case("$BOX_DIMENSIONS") !box dimensions
				do i=1,3
					read(this%fileUnit,*, iostat = ierr) this%box_dimensions(i)
					if(ierr .ne. 0) then
						print"(a,i1,a)", "BOXDATA ERROR: $BOX_DIMENSIONS(", i, ") does not contain proper data"
						error = .true.
					end if
				end do
				
			case("$LAYERS_LIMITS") !Z coordinates of layers
				!read up to 7 numbers...    
				do i=1,7
					read(this%fileUnit,*, iostat = ierr) line
					if(ierr .ne. 0) then
						print*, "problem when reading BOXDATA $LAYERS_LIMITS is probably the last line of the file"
					end if
					read(line,*, iostat = ierr) this%layers_limits(i)
					if(ierr .eq. 0) then
						this%layersCount = this%layersCount + 1
					else
						!no more numbers to read
						backspace(this%fileUnit, iostat = ierr)
						if(ierr .ne. 0) then
							print*, "backspace problem when reading $LAYERS_LIMITS"
							stop
						end if
						exit
					end if
				end do
						
			case("$INTERFACE_VOLUME_ELEMENT") !box dimensions
				do i=1,3
					read(this%fileUnit,*, iostat = ierr) this%interface_volume_element(i)
					if(ierr .ne. 0) then
						print"(a,i1,a)", "BOXDATA ERROR: $INTERFACE_VOLUME_ELEMENT(", i, ") does not contain proper data"
						error = .true.
					end if
				end do
				
			case("$FREQ") !frequency range for the fourier transform
				read(this%fileUnit,*, iostat = ierr) this%FREQ
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $FREQ does not contain proper data"
					error = .true.
				end if
				
			case("$DFREQ") !frequency range for the fourier transform
				read(this%fileUnit,*) this%DFREQ
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $NSTEP does not contain proper data"
					error = .true.
				end if

			case("$DT") !timestep in fs
				read(this%fileUnit,*, iostat = ierr) this%DT
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $DT does not contain proper data"
					error = .true.
				end if

			case("$CORRLEN") !length of correlation function in ps
				read(this%fileUnit,*, iostat = ierr) this%CORRLEN
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $CORRLEN does not contain proper data"
					error = .true.
				end if
			   
			case("$FILTER") !number of skipped frames while calculating correlation function
				if(this%filter == 0) then
					read(this%fileUnit,*, iostat = ierr) this%FILTER
					if(ierr .ne. 0) then
						print"(a)", "BOXDATA WARNING: $FILTER does not contain proper data"
						print*, "Filter will be set to default value"
						print*, ""
					end if
				else
					print*, "REMINDER:"
					print*, "$FILTER was modified by the program option"
					print*, "$FILTER value is set to:", this%filter
				end if
			
			case("$TEMPERATURE")
				read(this%fileUnit,*, iostat = ierr) this%TEMPERATURE
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $TEMPERATURE does not contain proper data"
					error = .true.
				end if

			case("$HBHIST_ANGLE_START")
				read(this%fileUnit,*, iostat = ierr) this%hbhist_anglestart
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $HBHIST_ANGLE_START does not contain proper data"
					error = .true.
				end if
				
			case("$HBHIST_ANGLE_END")
				read(this%fileUnit,*, iostat = ierr) this%hbhist_angleend
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $HBHIST_ANGLE_END does not contain proper data"
					error = .true.
				end if
				
			case("$HBHIST_DIST_START")
				read(this%fileUnit,*, iostat = ierr) this%hbhist_diststart
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $HBHIST_DIST_START does not contain proper data"
					error = .true.
				end if
				
			case("$HBHIST_DIST_END")
				read(this%fileUnit,*, iostat = ierr) this%hbhist_distend
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $HBHIST_DIST_END does not contain proper data"
					error = .true.
				end if
				
			case("$DIPOLE_R_START")
				read(this%fileUnit,*, iostat = ierr) this%dipole_rstart
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $DIPOLE_R_START does not contain proper data"
					error = .true.
				end if
				
			case("$DIPOLE_R_END")
				read(this%fileUnit,*, iostat = ierr) this%dipole_rend
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $DIPOLE_R_END does not contain proper data"
					error = .true.
				end if
				
			case("$DIPOLE_R_BINWIDTH")
				read(this%fileUnit,*, iostat = ierr) this%dipole_rbinwidth
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $DIPOLE_R_BINWIDTH does not contain proper data"
					error = .true.
				end if                
				
		   case("$DENSITY_MIN_R")
				read(this%fileUnit,*, iostat = ierr) this%density_min_r
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $DENSITY_MIN_R does not contain proper data"
					error = .true.
				end if
				
			case("$DENSITY_MAX_R")
				read(this%fileUnit,*, iostat = ierr) this%density_max_r
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $DENSITY_MAX_R does not contain proper data"
					error = .true.
				end if
				
			case("$DENSITY_TOL")
				read(this%fileUnit,*, iostat = ierr) this%density_tol
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $DENSITY_TOL does not contain proper data"
					error = .true.
				end if
						
			case("$NATOM") !number of atoms
				read(this%fileUnit,*, iostat = ierr) this%NATOM
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $NATOM does not contain proper data"
					error = .true.
				end if
					
			case("$NO") !number of oxygens
				read(this%fileUnit,*, iostat = ierr) this%NO
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $NATOM does not contain proper data"
					error = .true.
				end if
				
			case("$NSTEP") !length of trajectory
				read(this%fileUnit,*, iostat = ierr) this%NSTEP
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $NSTEP does not contain proper data"
					error = .true.
				end if
				
			case("$INTERFACE_SKIP") !number of skipped frames to calculate interface
				read(this%fileUnit,*, iostat = ierr) this%INTERFACE_SKIP
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $NSKIP does not contain proper data"
					error = .true.
				end if
				
			case("$CROSS_SKIP") !number of skipped frames while calculating correlation function
				read(this%fileUnit,*, iostat = ierr) this%CROSS_SKIP
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $CROSS_SKIP does not contain proper data"
					error = .true.
				end if
					
			case("$SELF_SKIP") !number of skipped frames while calculating correlation function
				read(this%fileUnit,*, iostat = ierr) this%SELF_SKIP
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $SELF_SKIP does not contain proper data"
					error = .true.
				end if
					
			case("$HBHIST_SKIP") !number of skipped frames while calculating correlation function
				read(this%fileUnit,*, iostat = ierr) this%HBHIST_SKIP
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $HBHIST_SKIP does not contain proper data"
					error = .true.
				end if
					
			case("$INTERFACE_PUSHBACK") !number of skipped frames while calculating correlation function
				read(this%fileUnit,*, iostat = ierr) this%INTERFACE_PUSHBACK
				
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $INTERFACE_PUSHBACK does not contain proper data"
					error = .true.
				end if
			
			case("$INTERFACE_DENSITY") !number of skipped frames while calculating correlation function
				read(this%fileUnit,*, iostat = ierr) this%INTERFACE_DENSITY
	 
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $INTERFACE_DENSITY does not contain proper data"
					error = .true.
				end if  
				
			case("$HBHIST_ANGLE_DIV")
				read(this%fileUnit,*, iostat = ierr) this%hbhist_anglediv
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $HBHIST_ANGLE_DIV does not contain proper data"
					error = .true.
				end if
				
			case("$HBHIST_DIST_DIV")    
				read(this%fileUnit,*, iostat = ierr) this%hbhist_distdiv
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $HBHIST_DIST_DIV does not contain proper data"
					error = .true.
				end if

			case("$DENSITY_BIN")
				read(this%fileUnit,*, iostat = ierr) this%density_bin
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $DENSITY_BIN does not contain proper data"
					error = .true.
				end if

			case("$HYDROXYL_METAL") !frequency range for the fourier transform
				read(this%fileUnit,'(A)') this%hydroxyl_metal
				if(ierr .ne. 0) then
					print"(a)", "BOXDATA ERROR: $HYDROXYL_METAL does not contain proper data"
					error = .true.
				end if

			case("$POLARIZATION")
				read(this%fileUnit,'(A)') line
				line = trim(line)
				line = adjustl(line)
				
				do i = 1, len(line)
					if(line(i:i) >= 'a' .and. line(i:i) <= 'z') then
						line(i:i) = achar(iachar(line(i:i)) - 32)
					end if
				end do
				
				select case(line)
					case ('SSP')
						this%P=1
						this%Q=1
						this%R=3
						this%polarization = "SSP"	

					case ('PPP')
						this%P=3
						this%Q=3
						this%R=3
						this%polarization = "PPP"	
					
					case default
						print*, "BOXDATA ERROR: polarization ", line, " is not supported."
						error = .true.
				end select

			case("$CANCEL_INTERFACE_CALCULATION")
				this%cancel_interface_calculation = .true.
				
			case("")
				
			case default
				print*, "REMAINDER: skipping ", trim(line), " - not a valid parameter"
				
			end select
				
		read(this%fileUnit,'(A)', IOSTAT = read_status) line
		
	end do
	
	if(error) then !if error in reading the file
		print*, ""
		print"(a)", "BOXDATA file has ERRORS listed above, please fix them"
		stop
	else
		print*, "reading of the BOXDATA file - OK"
		print*, ""
	end if
	
	!checking the input... catching errors in the initial values
	
	do i=1,3
		if(this%box_dimensions(i) <= 0) then
			print"(a,i1,a)", "ERROR: $BOX_DIMENSIONS(", i, ") is NOT positive nonzero value"
			error = .true.
		end if
	end do
	
	do i=1,3
		if(this%layers_limits(i) == 0) then
			print"(a,i1,a)", "WARNING: $LAYERS_LIMITS(", i, ") is zero"
			print*, "This will prevent Binder program to start"
		end if
	end do
	
	do i=1,3
		if(this%interface_volume_element(i) <= 0) then
			print"(a,i1,a)", "ERROR: $INTERFACE_VOLUME_ELEMENT(", i, ") is NOT positive nonzero value"
			error = .true.
		end if
	end do    
	
	if(this%DFREQ >= this%FREQ) then !dfreq is too high
		print*, "ERROR: $DFREQ is higher or equal to $FREQ"
		error = .true.
	end if
	
	if(this%DT <= 0) then
		print*, "ERROR: $DT is not positive nonzero value"
		error = .true.
	end if
	
	if(this%corrlen == 0) then
		print("(A)"), "WARNING:"
		print*, "$CORRLEN was not declared in BOXDATA file"
		print*, "Length of correlation function is set to $NSTEP in BOXDATA"
		print*, "This might lead to long calculation time"
		print*, ""
	else if(int(1000*this%CORRLEN/this%dt) > this%nstep) then
		print("(A)"), "WARNING:"
		print*, "$CORRLEN was longer than length of the simulation"
		print*, "$CORRLEN is set to the length of the simulation (ps): ", this%dt*this%nstep/1000
		print*, "This might lead to long calculation time"
		print*, ""
	end if
	
	if(this%TEMPERATURE <= 0) then
		print*, "ERROR: $TEMPERATURE is not positive nonzero value"
		error = .true.
	end if
	
	if(this%hbhist_angleend <= this%hbhist_anglestart) then
		print*, "ERROR: $HBHIST_ANGLE_END is lower or equal to $HBHIST_ANGLE_START"
		error = .true.
	end if
	
	if(this%hbhist_distend <= this%hbhist_diststart) then
		print*, "ERROR: $HBHIST_DIST_END is lower or equal to $HBHIST_DIST_START"
		error = .true.
	end if
	
	if(this%dipole_rend <= this%dipole_rstart) then
		print*, "ERROR: $DIPOLE_R_END is lower or equal to $DIPOLE_R_START"
		error = .true.
	end if
	
	if(this%dipole_rbinwidth <= 0) then
		print*, "ERROR: $DIPOLE_R_BINWIDTH is not positive nonzero value"
		error = .true.
	end if
	
	if(this%density_max_r <= this%density_min_r) then
		print*, "ERROR: $DENSITY_MAX_R is lower or equal to $DENSITY_MIN_R"
		error = .true.
	end if
	
	if(this%density_tol < 0) then
		print("(A)"), "WARNING:"
		print*, "$DENSITY_TOL is negative number - absolute value will be used"
		this%density_tol = -this%density_tol
	end if
	
	if(this%density_tol .eq. 0) then
		!adjust the density tol to be 1/(2*density_bin)
		this%density_tol = 1.0/(2.0*this%density_bin) 
	end if
	
	if(this%NATOM <= 0) then
		print("(A)"), "ERROR:"
		print*, "$NATOM is not positive nonzero value"
		error = .true.
	end if    
	
	if(this%NO <= 0) then
		print("(A)"), "ERROR:"
		print*, "$NO is not positive nonzero value"
		error = .true.
	end if
	
	if(this%NSTEP <= 0) then
		print("(A)"), "ERROR:"
		print*, "$NSTEP is not positive nonzero value"
		error = .true.
	end if    
	
	if(this%INTERFACE_SKIP < 0) then
		print("(A)"), "WARNING:"
		print*, "$INTERFACE_SKIP is negative number - absolute value will be used"
		this%INTERFACE_SKIP = -this%INTERFACE_SKIP
	end if  

	if(this%INTERFACE_SKIP == 0) then
		this%INTERFACE_SKIP = 1
	end if  

	if(this%INTERFACE_SKIP > this%NSTEP) then !cannot skip all the frames...
		print("(A)"), "ERROR:"
		print*, "$INTERFACE_SKIP is higher than the number of frames ($NSTEP) in BOXDATA"
		error = .true.
	end if
	
	if(this%CROSS_SKIP < 0) then
		print("(A)"), "WARNING:"
		print*, "$CROSS_SKIP is negative number - absolute value will be used"
		this%CROSS_SKIP = -this%CROSS_SKIP
	end if  

	if(this%CROSS_SKIP == 0) then
		this%CROSS_SKIP = 1
	end if  
	
	if(this%SELF_SKIP < 0) then
		print("(A)"), "WARNING:"
		print*, "$SELF_SKIP is negative number - absolute value will be used"
		this%SELF_SKIP = -this%SELF_SKIP
	end if  

	if(this%SELF_SKIP == 0) then
		this%SELF_SKIP = 1
	end if  

	if(this%HBHIST_SKIP < 0) then
		print("(A)"), "WARNING:"
		print*, "$HBHIST_SKIP is negative number - absolute value will be used"
		this%HBHIST_SKIP = -this%HBHIST_SKIP
	end if  

	if(this%HBHIST_SKIP == 0) then
		this%HBHIST_SKIP = 1
	end if  

	if (this%INTERFACE_PUSHBACK <= 0) then
		print("(A)"), "WARNING:"
		print*, "$INTERFACE_PUSHBACK in BOXDATA file is set to 0 or negative number"
		print*, "setting $INTERFACE_PUSHBACK to max value"
		print*, "this may lead to longer calculation time"
		this%INTERFACE_PUSHBACK = 100000000
	end if    
	
	if(this%hbhist_anglediv <= 0) then
		print("(A)"), "ERROR:"
		print*, "$HBHIST_ANGLE_DIV is not positive nonzero value"
		error = .true.
	end if
	
	if(this%hbhist_distdiv <= 0) then
		print("(A)"), "ERROR:"
		print*, "$HBHIST_DIST_DIV is not positive nonzero value"
		error = .true.
	end if

	if(this%density_bin <= 0) then
		print("(A)"), "ERROR:"
		print*, "$DENSITY_BIN is not positive nonzero value"
		error = .true.
	end if
	
	if(this%filter <= 0) then
		!filter <= 0 makes no sense -> Simones default version...
		!default filter value by Simone after the change of filter
		print("(A)"), "WARNING:"
		print*, "filter was set to a negative or zero value: ", this%filter
		this%filter = 1.06
		print*, "the filter value was modified to value: ", this%filter, " ps"
	else
		!filter is set according to the value in boxdata (can be overriden by -filter switch in program switches evaluation) 
	end if

	if(error) then
		print*, ""
		print*, "BOXDATA file has some ERRORS in the input values, please fix them"
		stop
	else
		print*, "checking the BOXDATA values - OK"
		print*, ""
	end if
	
	print*, "BOXDATA reading finished"
	print*, ""
	
	call this%print_recap()
	print*, ""
	
	close(this%fileUnit)
	
end subroutine read_boxdata

!!!returns number of samples of correlation fucntion based on DT and CORRLEN (ps)
integer*8 function get_maxlag(this)
	class(boxdata) ::                                       this
	
	if(this%corrlen == 0) then
		get_maxlag = this%nstep - 1
	else if(nint(1000*this%CORRLEN/this%dt) >= this%nstep) then
		get_maxlag = this%nstep - 1
	else
		get_maxlag = nint(1000*this%CORRLEN/this%DT)
	end if

end function get_maxlag

integer*8 function get_layersCount(this)
	class(boxdata) ::                                       this
	
	get_layersCount = this%layersCount
end function get_layersCount

subroutine print_recap(this)
	class(boxdata) ::										this
	integer ::												i
	
	print*, "BOXDATA RECAP:"
	print*, ""

	!arrays
	print "(A40,5X,F10.3)", adjustl("BOX_DIMENSIONS"), this%box_dimensions(1)
	print "(A40,5X,F10.3)", adjustl(""), this%box_dimensions(2)
	print "(A40,5X,F10.3)", adjustl(""), this%box_dimensions(3)
	print "(A40,5X,F10.3)", adjustl("INTERFACE_VOLUME_ELEMENT"), this%interface_volume_element(1)
	print "(A40,5X,F10.3)", adjustl(""), this%interface_volume_element(1)
	print "(A40,5X,F10.3)", adjustl(""), this%interface_volume_element(1)
	!todo
	print "(A40,5X,F10.3)", adjustl("LAYERS_LIMITS"), this%layers_limits(1)
	do i = 2, this%layersCount-1
		print "(A40,5X,F10.3)", adjustl(""), this%layers_limits(i)
	end do

	!reals
	print "(A40,5X,F10.3)", adjustl("FREQ"), this%FREQ
	print "(A40,5X,F10.3)", adjustl("DFREQ"), this%DFREQ
	print "(A40,5X,F10.3)", adjustl("DT"), this%DT
	print "(A40,5X,F10.3)", adjustl("CORRLEN"), this%CORRLEN
	print "(A40,5X,F10.3)", adjustl("FILTER"), this%FILTER
	print "(A40,5X,F10.3)", adjustl("TEMPERATURE"), this%TEMPERATURE
	print "(A40,5X,F10.3)", adjustl("HBHIST_ANGLESTART"), this%hbhist_anglestart
	print "(A40,5X,F10.3)", adjustl("HBHIST_ANGLE_END"), this%hbhist_angleend
	print "(A40,5X,F10.3)", adjustl("HBHIST_DIST_START"), this%hbhist_diststart
	print "(A40,5X,F10.3)", adjustl("HBHIST_DIST_END"), this%hbhist_distend
	print "(A40,5X,F10.3)", adjustl("DIPOLE_R_START"), this%dipole_rstart
	print "(A40,5X,F10.3)", adjustl("DIPOLE_R_END"), this%dipole_rend
	print "(A40,5X,F10.3)", adjustl("DIPOLE_R_BINWIDTH"), this%dipole_rbinwidth
	print "(A40,5X,F10.3)", adjustl("DENSITY_MIN_R"), this%density_min_r
	print "(A40,5X,F10.3)", adjustl("DENSITY_MAX_R"), this%density_max_r
	print "(A40,5X,F10.3)", adjustl("DENSITY_TOL"), this%density_tol
	print "(A40,5X,F10.3)", adjustl("INTERFACE_DENSITY"), this%interface_density

	!integers
	print "(A40,5X,I10)", adjustl("NATOM"), this%NATOM
	print "(A40,5X,I10)", adjustl("NO"), this%NO
	print "(A40,5X,I10)", adjustl("NSTEP"), this%NSTEP
	print "(A40,5X,I10)", adjustl("INTERFACE_SKIP"), this%INTERFACE_SKIP
	print "(A40,5X,I10)", adjustl("CROSS_SKIP"), this%CROSS_SKIP
	print "(A40,5X,I10)", adjustl("SELF_SKIP"), this%SELF_SKIP
	print "(A40,5X,I10)", adjustl("INTERFACE_PUSHBACK"), this%INTERFACE_PUSHBACK
	print "(A40,5X,I10)", adjustl("HBHIST_ANGLE_DIV"), this%hbhist_anglediv
	print "(A40,5X,I10)", adjustl("HBHIST_DIST_DIV"), this%hbhist_distdiv
	print "(A40,5X,I10)", adjustl("HBHIST_SKIP"), this%HBHIST_SKIP
	print "(A40,5X,I10)", adjustl("DENSITY_BIN"), this%density_bin

	!chars
	print "(A40,5X,A10)", adjustl("HYDROXYL_METAL"), adjustr(this%hydroxyl_metal)
	print "(A40,5X,A10)", adjustl("POLARIZATION"), adjustr(this%polarization)

	!logicals
	print "(A40,5X,L10)", adjustl("CANCEL_INTERFACE_CALCULATION"), this%cancel_interface_calculation
	

end subroutine print_recap

end module bdata
	
