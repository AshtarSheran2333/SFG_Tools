!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>                                                                           <!
!>                      SSP_CORR_P.f90                                       <!
!>                                                                           <!
!>  This program calculates SFG spectra from the MD simulations.             <!
!>  For reference see __GITHUB_LINK__                                        <!
!>  Or article __ARTICLE_LINK__                                              <!
!>                                                                           <!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program SSP_SFG_SELF_CROSS_CORR

use bdata
use SFG_FRAMES
use SFG_TYPES
use switches

implicit none

! boxdata collection
type(boxdata) ::                                            bd

! declaring the frame reader pointer and allocatable frame reader instances
class(frame_reader), pointer ::                             fr
type(trr_frame_reader), allocatable, target ::              trr
type(gro_frame_reader), allocatable, target ::              gro
type(xyz_frame_reader), allocatable, target ::              xyz

! MZ (molecule, time)
! AXX (molecule, time)
real*8, dimension(:,:), allocatable ::                      MZ,&
															AXX

! dAdRz polarizability derivatives
real*8, dimension(3,3) ::                                   dAdRz

! oxygen_position (r, molecule, time) - buffer for oxygen positions used in cross correlation calculations
real*4, dimension(:,:,:), allocatable ::                    oxygen_position

! av* - arrays storing average values of the MZ and AXX (molecule)
! *_corr (timelag) - arrays containing the correlation function
! crosscorr[AM] - holds cross correlation functions with different switch function, their average should be equal to cross_corr
real*8, dimension(:), allocatable ::                        avMZ,&
															avAXX,&
															hydroxyl_corr,&
															cross_corr,&
															self_corr,&
															AMSUM_corr,&
															crosscorrM, crosscorrA,&
															MM, A

! dMdRz - dipole moment derivatives
real*8, dimension(3) ::                                     dMdRz

! crossnorm - counter array to normalize the cross correlation function
! selfnorm - counter array to normalize the cross correlation function
integer*8, dimension(:), allocatable ::                     crossnorm,&
															selfnorm

! layers of interest hold the binder values to be included into the correlation function
integer*1, dimension(:), allocatable ::                     layers_of_interest

! binder_in_time holds maxlag frames of binder 
logical, dimension(:,:), allocatable ::                     binder_in_time

! timer - saves time to estimate run length
real*4 ::                                                   timer_1 = 0,&
															timer_2 = 0,&
															t_read = 0,&
															t_calc = 0,&
															timer_start = 0,&
															timer_estimate = 0

! ccp, ccpsq - cross correlation cutoff and its square
! dr - holds norm of vectors (midstep)
! scalar - holds dot product of vectors (midstep) used to enumerate D matrices
! auxiliary_variable - buffer for calculations...
real*8 ::                                                   ccp = 4,&
															ccpsq,&
															dr,&
															scalar,&
															auxiliary_variable,&
															benchmarkTime = 0

! t iterator through steps
! i,j,k,m,f general iterators
! timelag - lag of the correlation function
! nt1, nt0 - iterators used while calculating the correlation function
! maxlag - lag count to be evaluated in the correlation function
integer*8 ::                                                t,&
															i,&
															j,&
															k,&
															m,&
															n,&
															timelag,&
															corrlen,&
															nt1,&
															nt0

! ierr - to catch runtime errors
! narg - program argument count
! f* - file unit for *
! clk_* - used for timing...
integer ::                                                  ierr,&
															narg,&
															f_selfterms,&
															f_binder,&
															f_corr,&
															f_crossterms,&
															f_AMSUMterms,&
															f_hydroxylterms,&
															clk_1,clk_2,clk_rate,clk_start,clk_estimate

! nocheck - flag to ommit input file check
! enable_nlist - enable neighorlist when calculating crossterms
! skip_this - auxiliary variable...
logical ::                                                  nocheck = .false.,&
															skip_this,&
															enable_nlist = .false.,&
															benchmark = .false.

! binderFile - file name (binder file)
! output - file name (user defined prefix for the output files)
character(len=256) ::                                       binderFile = "binder.bin",&
															output = ""

! constants used by the program
real*8, parameter ::                                        pi = dacos(-1.d0),&
															kb = 1.386658d-23,&
															c = 2.99792458d10,& ! cm/s
															debye_to_ea = 0.208,& ! Debye to eA
															electron_to_coulomb = 1.6 ! 1e to C

! variables bound to the NLIST
real*8 ::                                                   maxdisplacement = 0,&
															secondmaxdisplacement = 0,&
															displacement,&
															rskin,&
															rskinsq,&
															maxdisplacementmaxt

! just a calculation buffer used to hold distance between two atoms
real*8, dimension(3) ::                                     diff

! nlist - holds idices of all neighbors of each atom
! neighbors - holds count of neighbors for each atom
integer*4, allocatable, dimension(:,:) ::                   nlist
integer*4, allocatable, dimension(:) ::                     neighbors

! strategy flags to make the conditions in the main loop simpler...
integer*1 ::                                                strategy_flags = b'00000000'

! those are all the possible computation strategies...
! if file contains hydroxyls, strategy flags should be set to b'10000000'
! if the main loop encounters any other flag... it should stop...
integer*1, parameter ::                                     sf_nothing =        b'00000000',&
															sf_self =           b'00000001',&
															sf_cross =          b'00000010',&
															sf_AMSUM =          b'00000100',&
															sf_hydroxyls =      b'10000000'

character(len = 256) :: op

!------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>                                                                           <!
!> values from Remi Khatib article for : dM(x,y,z)/dRz and dA(x,y,z)/dRz     <!
!> for water https://doi.org/10.1038/srep24287                               <!
!> where x,y,z are coordinates in the oh_ref                                 <!
!>                                                                           <!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dMdRz(1)   =-0.1500000     ! = dM(x)/dRz
dMdRz(2)   =-0.d0          ! = dM(y)/dRz
dMdRz(3)   = 2.1000000     ! = dM(z)/dRz
!------------------------------------------------------------------------------
dAdRz(1,1) = 0.4000000     ! = dA(x,x)/dRz
dAdRz(2,2) = 0.5300000     ! = dA(y,y)/dRz
dAdRz(3,3) = 1.5600000     ! = dA(z,z)/dRz

dAdRz(1,2) = 0.d0          ! = dA(x,y)/dRz
dAdRz(2,1) = dAdRz(1,2)    ! = dA(y,x)/dRz

dAdRz(1,3) = 0.0200000     ! = dA(x,z)/dRz
dAdRz(3,1) = dAdRz(1,3)    ! = dA(z,x)/dRz

dAdRz(2,3) = 0.d0          ! = dA(y,z)/dRz
dAdRz(3,2) = dAdRz(2,3)    ! = dA(z,y)/dRz
!------------------------------------------------------------------------------
! read BOXDATA

!this is extremely ugly...
if(iargc() < 1) then
	print*, "The program needs to be called with some arguments, see help -h"
	stop
end if

! is help switch present?
do i = 1, iargc()
	call getarg(i, op)
	if(getflag(trim(op)) == 'h') then
		call printHelp()
		stop
	end if
end do

!boxdata must be read first because of some switches...

call bd%read_boxdata()
if(.not. bd%hydroxyl_metal == "") then 
	! ugly way how to set the readers module global variable... whatever
	call set_hydroxyl_metal_name(bd%hydroxyl_metal)
end if

! evaluate program switches
call evaluate_switches()

! check wheather calculation strategies were set
if(strategy_flags == sf_nothing) then 
	print*, "no calculation strategy was set. Calculation is done."
	stop
end if

! checking the trajectory and binder files if applicable
if(.not. nocheck) then
	call check_input_files()
end if

! print recap of the parameters
call recap()

! open necessary files
call open_files()

! set the correlation function length + 1 for the lag 0
corrlen = bd%get_maxlag() + 1
	
! allocate all variables needed
call allocate_memory()

! start timing
if(benchmark == .true.) then
	call SYSTEM_CLOCK(clk_start, clk_rate)
end if
call SYSTEM_CLOCK(clk_2, clk_rate)

! start of the main calculation loop
! or rather frame reading loop :)
do t = 0, bd%nstep - 1
	! set the pointer of current time
	nt1 = mod(t,corrlen)
	! read frame
	call fr%read_frame()
	! select select molecules based on binder.bin and -L options
	if(.not. fr%is_hydroxyl()) then
		call fill_binder(f_binder, nt1)
	end if
	! timing of the reading
	call SYSTEM_CLOCK(clk_1, clk_rate)
	t_read = (clk_1-clk_2)/real(clk_rate)
	timer_1 = timer_1 + t_read

	! fill the AXX and MZ
	call fill_A_M(nt1)

	! every 1000th step processed print step
	if(mod(t,1000) == 0) print*, "step: ", t

	! adding the correlation function calculations should be as easy as this
	! the correlation function calculation strategies must be INDEPENDENT of each other...
	if(strategy_flags < 0) then 
		call corr_hydroxyls()
	else
		if(iand(strategy_flags, sf_self) == sf_self) call corr_self()
		if(iand(strategy_flags, sf_cross) == sf_cross) call corr_cross(enable_nlist)
		if(iand(strategy_flags, sf_AMSUM) == sf_AMSUM) call corr_AMSUM()
	end if
	
	! timing of the calculation
	call SYSTEM_CLOCK(clk_2, clk_rate)
	t_calc = (clk_2-clk_1)/real(clk_rate)
	timer_2 = timer_2 + t_calc

	if( benchmark ) then
		!record the time to maxlag (time to process grows with each frame)
		if( t == bd%get_maxlag() ) then
			call SYSTEM_CLOCK(clk_estimate, clk_rate)
		end if

		!do another X steps (3000) or minimum possible steps, time them properly
		if( t > bd%get_maxlag() .and. t < min((bd%get_maxlag() + 3000), (bd%get_maxlag() + min(3000, bd%nstep-bd%get_maxlag())))) then
			benchmarkTime = benchmarkTime + (t_read + t_calc)
		!terminate the program if X steps after maxlag passed
		else if (t > bd%get_maxlag()) then
			benchmarkTime = benchmarkTime + (t_read + t_calc)
			benchmarkTime = ((bd%nstep - bd%get_maxlag())/real(min(3000, bd%nstep-bd%get_maxlag()),8))*benchmarkTime
			benchmarkTime = benchmarkTime + (clk_estimate - clk_start)/real(clk_rate)
			print*, "BENCHMARK_ESTIMATION:"
			print*, "BE: ", benchmarkTime, " (s)"
			print*, "BE: ", benchmarkTime / 60, " (m)"
			print*, "BE: ", benchmarkTime / 3600, " (h)"
			print*, "BE: ", benchmarkTime / (3600 * 24), " (days)"
			print*, "The program execution stops now."
			stop
		end if
	end if

	! each 1k steps print estimated time and clear the timer
	if(mod(t,1000) == 0 .and. t > 0) then
		if(t < bd%get_maxlag()) then
			print*, "calibrating..., estimated time to finish ", s_to_HHHMMSS(((bd%nstep-1-t)/1000d0) * (timer_1 + timer_2))
		else
			print*, "estimated time to finish ", s_to_HHHMMSS(((bd%nstep-1-t)/1000d0) * (timer_1 + timer_2))
		end if
		print*, "Read/Calc time: ", timer_1/timer_2
		timer_1 = 0
		timer_2 = 0
	end if

	call backup()
end do

print*, ""
print*, "finishing the spectrum(a)..."
! normalize the correlation functions, do the unit conversion, write the corr files, make the spectra, and deallocate memory
! 1/10 ... just a factor to get the final result in 1e-22 m^2/V
if(strategy_flags < 0) then 
! hydroxyls
	do i=0,bd%get_maxlag()
		hydroxyl_corr(i)=hydroxyl_corr(i)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*(bd%nstep-i))
		write(f_hydroxylterms,*), i*bd%dt, hydroxyl_corr(i), bd%nstep-i
	end do
	close(f_hydroxylterms)
	call make_spectrum(hydroxyl_corr(:), trim(output)//"-hydroxyl")
	deallocate(hydroxyl_corr)
else
	if(iand(strategy_flags, sf_self) == sf_self) then
	! selfterms
		do i=0,bd%get_maxlag()
			self_corr(i)=self_corr(i)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*max(1,selfnorm(i)))
			write(f_selfterms,*), i*bd%dt, self_corr(i), max(1,selfnorm(i))
		end do
		deallocate(selfnorm)
		close(f_selfterms)
		call make_spectrum(self_corr(:), trim(output)//"-self")
		deallocate(self_corr)
	end if
	if(iand(strategy_flags, sf_cross) == sf_cross) then
		do i=0,bd%get_maxlag()
			cross_corr(i)=cross_corr(i)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*max(1,crossnorm(i)))
			crosscorrA(i)=crosscorrA(i)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*max(1,crossnorm(i)))
			crosscorrM(i)=crosscorrM(i)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*max(1,crossnorm(i)))
			write(f_crossterms,*), i*bd%dt, cross_corr(i), crosscorrM(i), crosscorrA(i), max(1,crossnorm(i))
		end do
		deallocate(crossnorm)
		close(f_crossterms)
		call make_spectrum(cross_corr(:), trim(output)//"-cross")
		call make_spectrum(crosscorrM(:), trim(output)//"-crossM")
		call make_spectrum(crosscorrA(:), trim(output)//"-crossA")
		deallocate(cross_corr)
		deallocate(crosscorrM)
		deallocate(crosscorrA)
	end if
	if(iand(strategy_flags, sf_AMSUM) == sf_AMSUM) then
		do i=0,bd%get_maxlag()
			AMSUM_corr(i)=AMSUM_corr(i)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*(bd%nstep-i))
			write(f_AMSUMterms,*), i*bd%dt, AMSUM_corr(i), bd%nstep-i
		end do
		close(f_AMSUMterms)
		call make_spectrum(AMSUM_corr(:), trim(output)//"-AMSUM")
		deallocate(AMSUM_corr)
	end if
end if

contains

subroutine backup()
	integer :: i, j, f, f_backup, ierr
	real :: wavenumber
	complex(kind=8), allocatable, dimension(:,:) :: intermediate_spectrum
	real(kind=8), allocatable, dimension(:) :: corr_buffer
	! 10 times per the trajectory length make spectra from currently calculated correlations
	! only the shifted spectra
	! only re and im
	! for cross skip cross[AM]
	! header based on the selected calculation strategies
	if( mod(t+1, bd%nstep/10) == 0 ) then
		! number of frequency points
		f = bd%FREQ/bd%DFREQ
		! allocate buffer for the spectra
		allocate(intermediate_spectrum(4, 0:f))
		! make the buffer zero
		intermediate_spectrum = 0
		! allocate buffer for the normalized corr function
		allocate(corr_buffer(0:bd%get_maxlag()))

		! open the backup file
		if( (t+1)/(bd%nstep/10) == 1 ) then
			! open the file in write mode when doing the first backup
			open(newunit = f_backup, file=trim(output)//"-backup.dat", recl = 160, iostat = ierr)
			if(ierr .ne. 0) then
				print*, "error opening", trim(output)//"-backup.dat"
				return
			end if
		else
			! open the file in append mode
			open(newunit = f_backup, file=trim(output)//"-backup.dat", access = 'append', recl = 160, iostat = ierr)
			if(ierr .ne. 0) then
				print*, "error opening", trim(output)//"-backup.dat"
				return
			end if
		end if
		
		! calculate FT to the buffers
		
		! write the file
		
		write(f_backup,*) "# timestep: ", t+1, " time [ps]: ", (t+1)*bd%dt/1000
		do i = 0, f
			wavenumber = dble(i)*bd%DFREQ ! cm^-1
			write(f_backup,"(f10.2,a)", advance = 'no') wavenumber, " "
			! self
			if(iand(strategy_flags, sf_self) == sf_self) then 
				if (i .eq. 0) then
					! convert units and normalize the corr function
					do j = 0, bd%get_maxlag()
						corr_buffer(j) = self_corr(j)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*max(1,selfnorm(j)))
					end do
					! calculate the spectrum
					call calc_spectrum(corr_buffer, intermediate_spectrum(1,:))
					intermediate_spectrum(1,:) = intermediate_spectrum(1,:) - intermediate_spectrum(1,f)
				end if
				write(f_backup,"(f16.8,a,f16.8,a)", advance = 'no') real(intermediate_spectrum(1,i)), " ", imag(intermediate_spectrum(1,i)), " "
			else
				write(f_backup,"(2I3)", advance = 'no') 0, 0
			end if
			! cross
			if(iand(strategy_flags, sf_cross) == sf_cross) then 
				if (i .eq. 0) then
					! convert units and normalize the corr function
					do j = 0, bd%get_maxlag()
						corr_buffer(j) = cross_corr(j)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*max(1,crossnorm(j)))
					end do
					! calculate the spectrum
					call calc_spectrum(corr_buffer, intermediate_spectrum(2,:))
					intermediate_spectrum(2,:) = intermediate_spectrum(2,:) - intermediate_spectrum(2,f)
				end if
				write(f_backup,"(f16.8,a,f16.8,a)", advance = 'no') real(intermediate_spectrum(2,i)), " ", imag(intermediate_spectrum(2,i)), " "
			else
				write(f_backup,"(2I3)", advance = 'no') 0, 0
			end if
			! AMSUM
			if(iand(strategy_flags, sf_AMSUM) == sf_AMSUM) then
				if (i .eq. 0) then
					! convert units and normalize the corr function
					do j = 0, bd%get_maxlag()
						corr_buffer(j) = AMSUM_corr(j)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*(t+1-j))
					end do
					! calculate the spectrum
					call calc_spectrum(corr_buffer, intermediate_spectrum(3,:))
					intermediate_spectrum(3,:) = intermediate_spectrum(3,:) - intermediate_spectrum(3,f)
				end if
				write(f_backup,"(f16.8,a,f16.8,a)", advance = 'no') real(intermediate_spectrum(3,i)), " ", imag(intermediate_spectrum(3,i)), " "
			else
				write(f_backup,"(2I3)", advance = 'no') 0, 0
			end if
			! hydroxyls
			if(iand(strategy_flags, sf_hydroxyls) == sf_hydroxyls) then
				if (i .eq. 0) then
					! convert units and normalize the corr function
					do j = 0, bd%get_maxlag()
						corr_buffer(j) = hydroxyl_corr(j)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*(t+1-j))
					end do
					! calculate the spectrum
					call calc_spectrum(corr_buffer, intermediate_spectrum(4,:))
					intermediate_spectrum(4,:) = intermediate_spectrum(4,:) - intermediate_spectrum(4,f)
				end if
				write(f_backup,"(f16.8,a,f16.8)", advance = 'no') real(intermediate_spectrum(4,i)), " ", imag(intermediate_spectrum(4,i))
			else
				write(f_backup,"(2I3)", advance = 'no') 0, 0
			end if
			write(f_backup,*) ""
		end do
		! end the block with extra two empty lines
		write(f_backup,*) ""
		write(f_backup,*) ""

		print*, "BACKUP", (t+1)/(bd%nstep/10)
		
		! close the file    
		close(f_backup)
		! deallocate buffer for the spectra
		deallocate(intermediate_spectrum)
	end if
end subroutine backup
	
! todo this is ugly...
subroutine calc_spectrum(in_corr, out_spectrum)
	real(kind=8), intent(IN), dimension(:) :: in_corr
	complex(kind=8), intent(INOUT), dimension(:) :: out_spectrum
	real :: wavenumber, omega, beta, time
	complex(kind=8) :: eiomegat, eiomegat1
	integer(kind=8) :: freq, lag, nf
	
	out_spectrum = 0
	
	nf = bd%FREQ/bd%DFREQ
	beta = 1/(kb*bd%TEMPERATURE)

	do freq = 0, nf
		wavenumber = dble(freq)*bd%DFREQ ! cm^-1
		omega = 2.d0*pi*c*wavenumber ! s^-1
		
		! trapezoidal integration
		do lag = 1, min(t+1,bd%get_maxlag()+1)
			time=bd%DT*dble(lag-1)*1.d-15 ! s
			eiomegat = dcmplx(dcos(dble(omega*time)), dsin(dble(omega*time)))
			time=bd%DT*dble(lag)*1.d-15 ! s
			eiomegat1 = dcmplx(dcos(dble(omega*time)), dsin(dble(omega*time)))
			if(lag < min(t+1,bd%get_maxlag()+1)) then
				out_spectrum(freq+1) = out_spectrum(freq+1) + 0.5d0 * ( eiomegat * in_corr(lag) * new_filter(lag-1, bd%filter, bd%dt) + eiomegat1 * in_corr(lag+1) * new_filter(lag, bd%filter, bd%dt) ) * bd%DT
			else ! add one artificial point to be 0 for the trap rule...
				out_spectrum(freq+1) = out_spectrum(freq+1) + 0.5d0 * eiomegat * in_corr(lag) * new_filter(lag-1, bd%filter, bd%dt) * bd%DT
			end if
		end do
		
		out_spectrum(freq+1) = dcmplx(0,-1) * beta * out_spectrum(freq+1) / omega
	end do
	

end subroutine calc_spectrum
	
! make_spectrum(arg1, arg2)
! arg1 - correlation function array
! arg2 - prefix of the spectrum file
! makes Fourier transform of correlation function, multiplies it by 1/kbT, and saves the spectrum to the file
! todo this is ugly
subroutine make_spectrum(corr, name)
	real(8), dimension(:), intent(IN) ::                    corr
	
	character(len=*), intent(IN) ::                         name
	
	complex(8), allocatable, dimension(:) ::                spectrum
	
	real(8) ::                                              wavenumber
	
	integer(8) ::                                           f,&
															f_spectrum,&
															f_shiftedspectrum
	
	f = bd%FREQ/bd%DFREQ

	allocate(spectrum(f+1))
	spectrum = 0
	
	call calc_spectrum(corr, spectrum)
	
	print*,""
	print*, name//"-spectrum.dat"
	open(newunit = f_spectrum, file=name//"-spectrum.dat", recl=120, iostat = ierr)
	if(ierr .ne. 0) then
		print"(a,a,a)", "ERROR: unable to create file '", name//"-spectrum.dat", "'"
		stop
	end if
	write(f_spectrum,"(A,F8.3)"), "#FILTER PARAMETER: ", bd%filter ! in ps
	write(f_spectrum,"(A)") "#wavelength cm-1, re, im, abs, phase deg"
	print*, name//"-shiftedspectrum.dat"
	print*, name//"-shiftedspectrum.dat is shifted by:"
	print*, "(Re, Im)", -spectrum(f)
	print*, "which is value of the spectrum at frequency ", f, "cm^-1"
	open(newunit = f_shiftedspectrum, file=name//"-shiftedspectrum.dat", recl=120, iostat = ierr)
	if(ierr .ne. 0) then
		print"(a,a,a)", "ERROR: unable to create file '", name//"-shiftedspectrum.dat", "'"
		stop
	end if
	write(f_shiftedspectrum,"(A,F8.3)"), "#FILTER PARAMETER: ", bd%filter ! in ps
	write(f_shiftedspectrum,"(A)") "#wavelength cm-1, re, im, abs, phase deg"
	
	do i=1,f+1
		! according to Khatib equation (3) gets the second order susceptibility
		wavenumber=dble(i-1)*bd%DFREQ       ! cm-1
	
		write(f_spectrum,*) int(wavenumber),&
					real(spectrum(i)),&
					imag(spectrum(i)),&
					abs(spectrum(i)),&
					datan(imag(spectrum(i))/real(spectrum(i)))*180/pi
	 
		! the spectrum is shifted by the value at frequency $FREQ, sine there should be no vibrations
		write(f_shiftedspectrum,*) wavenumber,&
					real(spectrum(i)-spectrum(f+1)),&
					imag(spectrum(i)-spectrum(f+1)),&
					abs(spectrum(i)-spectrum(f+1)),&
					datan(imag(spectrum(i)-spectrum(f+1))/real(spectrum(i)-spectrum(f+1)))*180/pi
	enddo
	
	close(f_spectrum)
	close(f_shiftedspectrum)
	
	deallocate(spectrum)
end subroutine make_spectrum

real*8 function new_filter(int_time, real_filter, dt)
	integer(kind=8), intent(IN) :: int_time ! steps
	real(kind=8), intent(IN) :: dt, real_filter !fs, ps

	new_filter = dexp(-(int_time*dt/(1000*real_filter))**2)

end function new_filter

real*8 function filter(time,tau)
	integer*8, intent(IN) ::                                time
	real*8, intent(IN) ::                                   tau
	
	filter = dexp(-(dble(time)/dble(tau))**2)
end
   
logical function in_range(molA, molB, timeA, timeB, cross_corr_parameter, cross_corr_parameter2)
	integer*8, intent(IN) :: molA, molB, timeA, timeB
	real*8, intent(IN) :: cross_corr_parameter, cross_corr_parameter2
	real*8, dimension(3) :: posDiff
	real*8 :: distance
	integer*1 :: i
	
	in_range = .FALSE.
	! calculate distance between two oxygens
	posDiff = oxygen_position(:,molA,timeA) - oxygen_position(:,molB,timeB)
	
	! PBC check
	do i=1,3
		if (posDiff(i) > bd%box_dimensions(i)/2) then
			posDiff(i) = posDiff(i) - bd%box_dimensions(i)
		else if (posDiff(i) < -bd%box_dimensions(i)/2) then
			posDiff(i) = bd%box_dimensions(i) + posDiff(i)
		end if
		if(posDiff(i) > cross_corr_parameter) return
	end do
	
	distance = posDiff(1)**2 + posDiff(2)**2 + posDiff(3)**2
	
	if(distance <= cross_corr_parameter2) in_range = .TRUE.
end function in_range

! evaluate_switches
! goes through program options and sets the program logic
subroutine evaluate_switches
	character(len=256) :: op, arg, arg1
	character :: option
	integer(kind=8) :: i,j
	real(kind=8) :: rl
	integer*1 ::  layerBoundary_1 = 100,&
				layerBoundary_2 = 100,& 
				swap_i1 = 100
	logical :: OK = .false.

	if(iargc() < 1) then
		print*, "The program needs to be called with some arguments, see help -h"
		stop
	end if

	! is help switch present?
	do i = 1, iargc()
		call getarg(i, op)
		if(getflag(trim(op)) == 'h') then
			call printHelp()
			stop
		end if
	end do

	i = 1
	do
		call getarg(i, op)
		option = getFlag(trim(op))

		select case(option)
			case('i')
				call getSwitchString(i, op, arg)
				! INPUT GRO
				if(index(arg, '.gro') .ne. 0) then
					print*, "input file:"
					print*, ""
					print"(tr5,a,tr5,a)", ".gro multiple frame file: ", trim(arg)
					print*, ""
					allocate(gro)
					fr => gro
					call fr%open_file(arg)
				! INPUT XYZ
				else if(index(arg, '.xyz') .ne. 0) then
					call getSwitchString(i, op, arg1)
					if(index(arg1, '.xyz') .ne. 0) then
						print*, "input files:"
						print*, ""
						print"(tr5,a,tr5,a)", ".xyz position file: ", trim(arg)
						print"(tr5,a,tr5,a)", ".xyz velocity file: ", trim(arg1)
						print*, ""
						allocate(xyz)
						fr => xyz
						call fr%open_file(arg, arg1)
					else
						print*, "program needs position and velocity .xyz files in order position file velocity file"
						stop
					end if
				! INPUT TRR
				else if(index(arg, '.trr') .ne. 0) then
					call getSwitchString(i, op, arg1)
					if(index(arg1, '.gro') .ne. 0) then
						print*, "input files:"
						print*, ""
						print"(tr5,a,tr5,a)", ".trr trajectory file: ", trim(arg)
						print"(tr5,a,tr5,a)", ".gro atleast single frame file: ", trim(arg1)
						print*, ""
						allocate(trr)
						fr => trr
						call fr%open_file(arg, arg1)
					else
						print*, "program needs .trr and .gro files in order .trr file .gro file"
						stop
					end if
				else
					print*, "Unknown input file format: ", trim(arg)
					stop
				end if

				if(fr%is_hydroxyl()) then
					strategy_flags = ior(strategy_flags, sf_hydroxyls)
				end if
				i = i + 1
			case('l')
				call getSwitchInteger(i, op, layerBoundary_1)
				call getSwitchInteger(i, op, layerBoundary_2)

				!layers can be: TOP | 8 ... 15 | 7 ... 0 | BOTTOM
				!so basically combinations like 0 7, 1 3 will do only that interval
				!likewise 8 15, 9 13
				!problem is when both interfaces are specified:
				! 0 15 should do 0 ... 7 15
				! 7 9 should do 7 9 ... 15
				! the layer boundaries needs to be sorted
				
				!sort layer boundaries
				if(layerBoundary_1 > layerBoundary_2) then
					swap_i1 = layerBoundary_1
					layerBoundary_1 = layerBoundary_2
					layerBoundary_1 = swap_i1
				end if
				!both boundaries are < 7 ? do interval
				!both boundaries are > 7 ? do interval
				if( ((layerBoundary_1 < 8) .and. (layerBoundary_2 < 8)) .or. ((layerBoundary_1 > 7) .and. (layerBoundary_2 > 7)) ) then
					allocate(layers_of_interest(layerBoundary_2 - layerBoundary_1 + 1))
					do j = layerBoundary_1, layerBoundary_2
						layers_of_interest(j+1-layerBoundary_1) = j
					end do
				!otherwise, list lower boundary to 7, and upper boundary to 15
				else
					j = (8 - layerBoundary_1) + (16 - layerBoundary_2)
					allocate(layers_of_interest(j))
					do j = layerBoundary_1, 7
						layers_of_interest(j+1-layerBoundary_1) = j
					end do
					do j = layerBoundary_2, 15
						layers_of_interest( (j+1-layerBoundary_2) + (8 - layerBoundary_1) ) = j
					end do
				end if
				
				print*, "selected layer binder flags: ", layers_of_interest(:)
				print*, ""
				i = i + 1
			case('o')
				call getSwitchString(i, op, output)
				print*, "output file prefix: ", trim(output)
				print*, ""
				i = i + 1
			case('s')
				strategy_flags = ior(strategy_flags, sf_self)
				i = i + 1
			case('c')
				strategy_flags = ior(strategy_flags, sf_cross)
				call getSwitchReal(i, op, ccp)
				ccpsq = ccp * ccp
				i = i + 1
			case('n')
				enable_nlist = .true.
				i = i + 1
			case('a')
				strategy_flags = ior(strategy_flags, sf_AMSUM)
				i = i + 1
			case('d')
				call getSwitchString(i, op, arg)
				open(99, FILE=trim(arg), iostat=ierr, status = "old")
				
				if(ierr .ne. 0) then
					print"(a,a,a)", "ERROR: -D file '", trim(arg), "' does not exist"
					stop
				end if
				
				do j=1,9
					read(99,*, iostat = ierr) rl
					if(ierr .ne. 0) then
						print"(a,a,a,i1,a)", "ERROR: -D ", trim(arg), " line ", j, " is in incorrect format"
						stop
					end if
				end do
				
				rewind(99)
				
				read(99,*) dMdRz(1) ! = dM(x)/dRz
				read(99,*) dMdRz(2) ! = dM(y)/dRz
				read(99,*) dMdRz(3) ! = dM(z)/dRz
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				read(99,*) dAdRz(1,1) ! = dA(x,x)/dRz
				read(99,*) dAdRz(2,2) ! = dA(y,y)/dRz
				read(99,*) dAdRz(3,3) ! = dA(z,z)/dRz
				
				read(99,*) dAdRz(1,2) ! = dA(x,y)/dRz
				dAdRz(2,1) = dAdRz(1,2) ! = dA(y,x)/dRz
				
				read(99,*) dAdRz(1,3) ! = dA(x,z)/dRz
				dAdRz(3,1) = dAdRz(1,3) ! = dA(z,x)/dRz
				
				read(99,*) dAdRz(2,3) ! = dA(y,z)/dRz
				dAdRz(3,2) = dAdRz(2,3) ! = dA(z,y)/dRz

				close(99)
				i = i + 1
			case('g')
				nocheck = .true.
				i = i + 1
			case('f')
				call getSwitchReal(i, op, rl)
				if(rl > 0) then
					print*, "BOXDATA $FILTER value was overriden by -F(filter) [ps]", rl
					bd%filter = rl
				else
					print*, "WARNING: -F(filter) <arg> has a invalid value, BOXDATA $FILTER will be applied [ps]: ", bd%filter
				end if
				i = i + 1
			case('b')
				benchmark = .true.
				i = i + 1
			case default
				print*, "skipping invalid option ", trim(op)
				i = i + 1
		end select

		! all the switches were evaluated
		if(i > iargc()) exit
	end do

	! check if all the mandatory options were selected
	if(.not. associated(fr)) then
		print*, "input file(s) must be specified, see help -h"
		stop
	end if
	! layers does not have to be specified when dealing with hydroxyls
	if(.not. (allocated(layers_of_interest) .or. fr%is_hydroxyl())) then
		print*, "the layers to be analysed must be specified, see help -h"
		stop
	end if

	if(fr%is_hydroxyl()) then
		strategy_flags = sf_hydroxyls
	end if

end subroutine evaluate_switches

! printHelp()
! prints help dialog of the program
subroutine printHelp()
	! this is the help dialog...

	print "(a)",        ""
	print "(a)",        "__________________________Help dialog of SSP_CORR______________________________"
	print "(a)",        ""                          
	print "(a)",        "_____________________________Mandatory options_________________________________"
	print "(a)",        ""
	print "(a)",        "-I(input) <file1> *<file2>"
	print "(tr4,a)",        "specifies the input files. Options are:"
	print "(tr8,a)",            "*.xyz (positions), *.xyz (velocities)"
	print "(tr8,a)",            "*.gro (positions+velocities)"
	print "(tr8,a)",            "*.trr, *.gro (one frame - can be only positions)"
	print "(a)", ""
	print "(a)",        "-L(layer) <number1 (from)> <number2 (to)>"
	print "(tr4,a)",        "specifies the layers to be analyzed. Options are:"
	print "(tr8,a)",            "BOTTOM | 0 1 2 3 4 5 6 7 | 15 14 13 12 11 10 9 8 | TOP"
	print "(a)",        ""
	print "(a)",        "___________________________Optional useful options_____________________________"
	print "(a)",        ""
	print "(a)",        "-O(output) <name>"
	print "(tr4,a)",        "<name> will be appended by the output file names"
	print "(a)",        ""
	print "(a)",        "-S(self)"
	print "(tr4,a)",        "enables calculation of the cross correlation terms"
	print "(a)",        ""
	print "(a)",        "-C(cross) <real number>"
	print "(tr4,a)",        "enables calculation of the cross correlation terms"
	print "(tr4,a)",        "<real number> is the cutoff radius for the cross correlation terms in [A]" 
	print "(a)",        ""
	print "(a)",        "-N(neighborlist)"
	print "(tr4,a)",        "enables neighbor list for for the cross correlation calculation"
	print "(a)",        ""
	print "(a)",        "-A(AMSUM)"
	print "(tr4,a)",        "enables calculation of the correlation function of summed dipole moments and polarizabilities"
	print "(a)",        ""
	print "(a)",        "-D(dipole_moment_and_polarizability_derivatives) <filename>"
	print "(tr4,a)",        "file with dipole moment and polarizability derivatives"
	print "(tr4,a)",        "file format is described in the documentation"
	print "(tr4,a)",        "default values are described in the documentation"
	print "(a)",        ""
	print "(a)",        "-F(filter) <real number>"
	print "(tr4,a)",        "filter parameter in [ps]"
	print "(tr4,a)",        "overrides $FILTER in BOXDATA file"
	print "(tr4,a)",        "filter function is described in the documentation"
	print "(a)",        ""
	print "(a)",        "-G(good)"
	print "(tr4,a)",        "skips checking of the input files"
	print "(a)",        ""
	print "(a)",        "-B(benchmark)"
	print "(tr4,a)",        "runs a benchmark (processes corrlen steps + up to 3000 steps)"
	print "(tr4,a)",		"then estimates the time to finish the calculation"
	print "(a)",        ""
	print "(a)",        "-H(help)"
	print "(tr4,a)",        "prints this help dialog"
	print "(a)",        ""
	print "(a)",        "for more details see documentation"
end subroutine printHelp
	
! cross_product(arg1, arg2)
! arg1, arg2 - 3d array vector
! returns cross product of arg1 and arg2
function cross_product(a,b)
	real*8, dimension(3) :: cross_product
	real*8, dimension(3), intent(IN) :: a, b

	cross_product(1) = a(2) * b(3) - a(3) * b(2)
	cross_product(2) = a(3) * b(1) - a(1) * b(3)
	cross_product(3) = a(1) * b(2) - a(2) * b(1)
end function cross_product

subroutine get_binder_frame(channel, binder_frame)
	integer, intent(IN) :: channel
	integer*1, dimension(:), intent(OUT) :: binder_frame
	do i=1,bd%NO,2
		read(channel, iostat = ierr) binder_frame(i)
		if(ierr .ne. 0) then
			print*, "ERROR: 'binder.bin' does not have enough data inside"
			stop
		end if
		if((i+1) <= bd%no) then
			binder_frame(i+1) = binder_frame(i) .AND. 15 ! lower nibble
		end if
		binder_frame(i) = binder_frame(i) / 16 ! upper nibble
	end do
end subroutine get_binder_frame

function fit_in_pbc(vector, bd)
	real*8, dimension(3) :: fit_in_pbc, vector, bd
	integer :: i
	
	do i=1,3
		do while (vector(i) > bd(i)/2)
			vector(i) = vector(i) - bd(i)
		end do
		
		do while (vector(i) < -bd(i)/2)
			vector(i) = bd(i) + vector(i)
		end do
	end do
	
	fit_in_pbc = vector
end function fit_in_pbc

subroutine recap()
	print*, ""
	print*, "Recap of calculation settings:"
	print*, ""
	if(output .ne. "") then
		print*, "Output file prefix: ", trim(output)
		print*, ""
	else
		print*, "Output file prefix: ", "-no prefix selected-"
		print*, ""
	end if
	
	if(iand(strategy_flags, sf_self) .eq. sf_self) then ! self
		print*, "self correlation terms will be calculated"
		print*, "self correlation function"
		print*, "will be stored in file: <output file prefix>-selfterms.dat"
		print*, "self correlation spectrum"
		print*, "will be stored in files: <output file prefix>-self-*.dat"
		print*, ""
	end if

	if(iand(strategy_flags, sf_cross) .eq. sf_cross) then ! cross
		print*, "cross correlation terms will be calculated"
		if(enable_nlist .eqv. .true.) then
			print*, "neighborlist will be used"
		else
			print*, "neighbor list will not be used"
		end if
		print*, "r_cut: ", ccp, "Angstrom"
		print*, "cross correlation function"
		print*, "will be stored in file: <output file prefix>-crossterms.dat"
		print*, "cross correlation spectrum"
		print*, "will be stored in files: <output file prefix>-cross-*.dat"
		print*, ""
	end if

	if(iand(strategy_flags, sf_AMSUM) .eq. sf_AMSUM) then ! AMSUM
		print*, "corelation of whole selected layer summed properties will be calculated"
		print*, "the correlation function"
		print*, "will be stored in file: <output file prefix>-AMSUMterms.dat"
		print*, "the spectrum"
		print*, "will be stored in file: <output file prefix>-AMSUMterms.dat"
		print*, ""
	end if

	print*, "dt is set to [fs]: ", bd%DT
	print*, "correlation function length [ps]:", bd%get_maxlag()*bd%dt/1000
	print*, "temperature is set to [K]: ", bd%TEMPERATURE
	print*, "Filter parameter is set to [ps]: ", bd%filter ! ps
	print*, "" ! dM/dr dA/dr
	print*, "dM_x/dr_z: ", dMdRz(1)
	print*, "dM_x/dr_z: ", dMdRz(2)
	print*, "dM_x/dr_z: ", dMdRz(3)
	print*, ""
	print*, "dA_xx/dr_z: ", dAdRz(1,1)
	print*, "dA_yy/dr_z: ", dAdRz(2,2)
	print*, "dA_zz/dr_z: ", dAdRz(3,3)
	print*, "dA_xy/dr_z: ", dAdRz(1,2)
	print*, "dA_xz/dr_z: ", dAdRz(1,3)
	print*, "dA_yz/dr_z: ", dAdRz(2,3)
	
	print*, "_______________________________________________________________________________"
end subroutine recap

subroutine check_input_files()
	integer*1, dimension(:), allocatable ::                     binder
	allocate(binder(bd%NO))
	
	! CAREFUL ... working with files is clumsy that way...
	print*, "Checking the contents of trajectory file(s)"
	print*, ""
	do i = 1, bd%NSTEP
		call fr%skip_frame()
	end do
	call fr%rewind_file()
	print*, "contents of trajectory files - OK"
	print*, ""

	if(.not. fr%is_hydroxyl()) then
		print*, "Checking the 'binder.bin' file"
	
		open(newunit = f_binder, file = binderFile, form = "unformatted", access = 'stream', convert = 'big_endian', status = 'old', iostat = ierr)
		if(ierr .ne. 0) then
			print"(a,a,a)", "ERROR: unable to open file '", binderFile, "'"
			stop
		end if
	
		do t=1, bd%NSTEP
			call get_binder_frame(f_binder, binder)
		end do
	
		close(f_binder)
	
		print*, "'binder.bin' file - OK"
		print*, ""
	end if
	
	deallocate(binder)
	print*, "_______________________________________________________________________________"
end subroutine check_input_files

subroutine open_files()
	! todo this is ugly, but someone is too lazy to wrap this...

	if(strategy_flags < 0) then 
		! initialize -hydroxylterms.dat
		open(newunit = f_hydroxylterms, FILE=trim(output)//'-hydroxylterms.dat', iostat = ierr)
		if(ierr .ne. 0) then
			print"(a,a,a)", "ERROR: unable to create file '", trim(output)//'-hydroxylterms.dat',"'"
			stop
		end if
		! file header
		write(f_hydroxylterms,*) "#N frames: ", bd%nstep
		write(f_hydroxylterms,*) "#DT: ", bd%dt
		write(f_hydroxylterms,*) "#layers: none - hydroxyls"
		write(f_hydroxylterms,"(A,9F8.4)") "#M and A derivatives: ", dMdRz(1),& ! = dM(x)/dRz
			dMdRz(2),& ! = dM(y)/dRz
			dMdRz(3),& ! = dM(z)/dRz
			dAdRz(1,1),&  ! = dA(x,x)/dRz
			dAdRz(2,2),&  ! = dA(y,y)/dRz
			dAdRz(3,3),&  ! = dA(z,z)/dRz
			dAdRz(1,2),&  ! = dA(x,y)/dRz
			dAdRz(1,3),&  ! = dA(x,z)/dRz
			dAdRz(2,3) ! = dA(y,z)/dRz
		write(f_hydroxylterms,"(A)") "#timelag, hydroxyls_corr_func, norm"  
	else
		! if not hydroxyls, binder is needed...
		open(newunit = f_binder, file = binderFile, form = "unformatted", access = 'stream', convert = 'big_endian', status = 'old', iostat = ierr)
		if(ierr .ne. 0) then
			print"(a,a,a)", "ERROR: unable to open file '", binderFile, "'"
			stop
		end if
		
		! strategy self
		if(iand(strategy_flags, sf_self) == sf_self) then
			! initialize the -selfterms.dat file
			open(newunit = f_selfterms, FILE=trim(output)//'-selfterms.dat', iostat = ierr)
			if(ierr .ne. 0) then
				print"(a,a,a)", "ERROR: unable to create file '", trim(output)//'-selfterms.dat',"'"
				stop
			end if
			! file header
			write(f_selfterms,*) "#N frames: ", bd%nstep
			write(f_selfterms,*) "#DT: ", bd%dt, " fs"
			write(f_selfterms,*) "#SELF_SKIP: ", bd%self_skip, " frames"
			write(f_selfterms,*) "#layers (binder flags): ", layers_of_interest(:)
			write(f_selfterms,"(A,9F8.4)") "#M and A derivatives: ", dMdRz(1),& ! = dM(x)/dRz
				dMdRz(2),& ! = dM(y)/dRz
				dMdRz(3),& ! = dM(z)/dRz
				dAdRz(1,1),&  ! = dA(x,x)/dRz
				dAdRz(2,2),&  ! = dA(y,y)/dRz
				dAdRz(3,3),&  ! = dA(z,z)/dRz
				dAdRz(1,2),&  ! = dA(x,y)/dRz
				dAdRz(1,3),&  ! = dA(x,z)/dRz
				dAdRz(2,3) ! = dA(y,z)/dRz
			write(f_selfterms,"(A)") "#timelag, corr_func, norm"  
		end if
		
		! strategy cross
		if(iand(strategy_flags, sf_cross) == sf_cross) then
			! initialize the -crossterms file
			open(newunit = f_crossterms, FILE=trim(output)//'-crossterms.dat', iostat = ierr, recl=160)
			if(ierr .ne. 0) then
				print"(a,a,a)", "ERROR: unable to create file '", trim(output)//'-crossterms.dat',"'"
				stop
			end if
			! file header
			write(f_crossterms,*) "#N frames: ", bd%nstep
			write(f_crossterms,*) "#DT: ", bd%dt, " fs"
			write(f_crossterms,*) "#CROSS SKIP: ", bd%cross_skip
			write(f_crossterms,*) "#Cross cutoff radius: ", ccp, " Angstrom"
			write(f_crossterms,*) "#layers (binder flags): ", layers_of_interest(:)
			write(f_crossterms,"(A,9F8.4)") "#M and A derivatives: ", dMdRz(1),& ! = dM(x)/dRz
				dMdRz(2),& ! = dM(y)/dRz
				dMdRz(3),& ! = dM(z)/dRz
				dAdRz(1,1),&  ! = dA(x,x)/dRz
				dAdRz(2,2),&  ! = dA(y,y)/dRz
				dAdRz(3,3),&  ! = dA(z,z)/dRz
				dAdRz(1,2),&  ! = dA(x,y)/dRz
				dAdRz(1,3),&  ! = dA(x,z)/dRz
				dAdRz(2,3) ! = dA(y,z)/dRz
			write(f_crossterms,"(A)") "#timelag, cross_corr_func, cross_corr_funcM, cross_corr_funcA, norm"  
		end if
		
		! strategy AMSUM
		if(iand(strategy_flags, sf_AMSUM) == sf_AMSUM) then
			! initialize -AMSUMterms.dat
			open(newunit = f_AMSUMterms, FILE=trim(output)//'-AMSUMterms.dat', iostat = ierr)
			if(ierr .ne. 0) then
				print"(a,a,a)", "ERROR: unable to create file '", trim(output)//'-hydroxylterms.dat',"'"
				stop
			end if
			! file header
			write(f_AMSUMterms,*) "#N frames: ", bd%nstep
			write(f_AMSUMterms,*) "#DT: ", bd%dt, " fs"
			write(f_AMSUMterms,*) "#layers (binder flags): ", layers_of_interest(:)
			write(f_AMSUMterms,"(A,9F8.4)") "#M and A derivatives: ", dMdRz(1),& ! = dM(x)/dRz
				dMdRz(2),& ! = dM(y)/dRz
				dMdRz(3),& ! = dM(z)/dRz
				dAdRz(1,1),&  ! = dA(x,x)/dRz
				dAdRz(2,2),&  ! = dA(y,y)/dRz
				dAdRz(3,3),&  ! = dA(z,z)/dRz
				dAdRz(1,2),&  ! = dA(x,y)/dRz
				dAdRz(1,3),&  ! = dA(x,z)/dRz
				dAdRz(2,3) ! = dA(y,z)/dRz
			write(f_AMSUMterms,"(A)") "#timelag, AMSUM_corr_func, norm"  
		end if
	end if
end subroutine open_files

subroutine allocate_memory()
	! allocate and set general variables
	allocate(Axx(bd%NO, 0:bd%get_maxlag()))
	allocate(Mz(bd%NO, 0:bd%get_maxlag()))
	allocate(avMz(bd%NO))
	allocate(avAXX(bd%NO))
	axx = 0
	mz = 0
	avMz=0
	avAXX=0

	! strategy hydroxyls
	if(strategy_flags < 0) then 
		allocate(hydroxyl_corr(0:bd%get_maxlag()))
		hydroxyl_corr = 0
	else
		! if not dealing with hydroxyls, binder is needed
		allocate(binder_in_time(BD%NO,0:bd%get_maxlag()))

		! strategy self
		if(iand(strategy_flags, sf_self) == sf_self) then
			allocate(self_corr(0:bd%get_maxlag()))
			allocate(selfnorm(0:bd%get_maxlag()))
			selfnorm = 0
			self_corr = 0
		end if
		
		! strategy cross
		if(iand(strategy_flags, sf_cross) == sf_cross) then
			allocate(cross_corr(0:bd%get_maxlag()))
			allocate(crosscorrM(0:bd%get_maxlag()))
			allocate(crosscorrA(0:bd%get_maxlag()))
			crosscorrM = 0
			crosscorrA = 0
			allocate(oxygen_position(3,bd%NO,0:bd%get_maxlag()))
			cross_corr = 0
			! NLIST
			allocate(nlist(bd%NO,bd%NO))
			nlist = 0
			! CAREFUL first dimension can be shrinked...
			! this did not cause any problems since our systems are < 10k water molecules
			! this could be done easily in c++ utilizing std::vector...
			! I dont want to implement dynamic array RN since this is just experimental feature
			! I can theoretically add another calculation strategy for cross_nlist... 
			! to save all this space... and some computation costs
			allocate(neighbors(bd%NO))
			neighbors = 0
			allocate(crossnorm(0:bd%get_maxlag()))
			crossnorm = 0
			
			if(.not. enable_nlist) then
				! just fill the neighbor list with all possible combinantions
				do m = 1, bd%NO
					do n = 1, bd%NO
						if(m .eq. n) cycle
						neighbors(m) = neighbors(m) + 1
						nlist(neighbors(m),m) = n
					end do
				end do
			end if
		end if
		
		! strategy AMSUM
		if(iand(strategy_flags, sf_AMSUM) == sf_AMSUM) then
			allocate(A(0:bd%get_maxlag()))
			A = 0
			allocate(MM(0:bd%get_maxlag()))
			MM = 0
			allocate(AMSUM_corr(0:bd%get_maxlag()))
			AMSUM_corr = 0
		end if
	end if
	
end subroutine allocate_memory

subroutine fill_A_M(time_index)
	integer*8, intent(IN) :: time_index
	integer :: m,n,i,j,b=2
	real*8 :: dr, scalar
	real*8, dimension(3) :: h1, h2
	real*8, dimension(2,3,3) :: D
	real*8, dimension(2) :: vz

	Axx(:,time_index) = 0
	Mz(:,time_index) = 0
	
	avMz=0.d0
	avAXX=0.d0

	! print*, "evaluating D matrices and v(z) for each mol"
	!$OMP PARALLEL DEFAULT(SHARED)
	!$OMP DO PRIVATE(m,n,i,j,h1,h2,dr,scalar,D,vz)
	do m=1,bd%NO
		
		!!!!!!!!!!! z(1) !!!!!!!!!!
		h1 = fr%molecule(m)%h1%position - fr%molecule(m)%o%position
		! pbc check
		do i=1,3
			do while (h1(i) > bd%box_dimensions(i)/2d0)
				h1(i) = h1(i) - bd%box_dimensions(i)
			end do
			
			do while (h1(i) < -bd%box_dimensions(i)/2d0)
				h1(i) = bd%box_dimensions(i) + h1(i)
			end do
		end do
		dr = norm2(h1)
		D(1,:,3) = h1/dr
		!!!!!!!!!!! z(2) !!!!!!!!!!
		h2 = fr%molecule(m)%h2%position - fr%molecule(m)%o%position
		! pbc check
		do i=1,3
			do while (h2(i) > bd%box_dimensions(i)/2d0)
				h2(i) = h2(i) - bd%box_dimensions(i)
			end do
			
			do while (h2(i) < -bd%box_dimensions(i)/2d0)
				h2(i) = bd%box_dimensions(i) + h2(i)
			end do
		end do
		dr = norm2(h2)
		D(2,:,3) = h2/dr
		!!!!!!!!!!! x(1) !!!!!!!!!!
		scalar = dot_product(h2,D(1,:,3))
		D(1,:,1) = scalar*D(1,:,3)-h2
		dr = norm2(D(1,:,1))
		if(dr == 0) then 
			if((D(1,1,3) == 0) .and. (D(1,2,3) == 0)) then
				! be careful, in this case we would divide by zero!
				! X axis if lab Z == Z else -X axis
				D(1,:,1) = (/D(1,3,3),0d0,0d0/)
			else
				! if hoh angle is 180deg (should be impossible in case of water molecules - but useful in case of surface hydroxyls)
				! chose X perpendicular to Z
				D(1,:,1) = (/D(1,2,3),-D(1,1,3),0d0/)
				dr = norm2(D(1,:,1))
				D(1,:,1) = D(1,:,1)/dr
			end if
		else
			D(1,:,1) = D(1,:,1)/dr
		end if
		!!!!!!!!!!! x(2) !!!!!!!!!!
		scalar = dot_product(h1,D(2,:,3))
		D(2,:,1) = scalar*D(2,:,3)-h1
		dr = norm2(D(2,:,1))
		if(dr == 0) then 
			if((D(2,1,3) == 0) .and. (D(2,2,3) == 0)) then
				! be careful, in this case we would divide by zero!
				! X axis if lab Z == Z else -X axis
				D(2,:,1) = (/D(2,3,3),0d0,0d0/)
			else
				! if hoh angle is 180deg (should be impossible in case of water molecules - but useful in case of surface hydroxyls)
				! chose X perpendicular to Z
				D(2,:,1) = (/D(2,2,3),-D(2,1,3),0d0/)
				dr = norm2(D(2,:,1))
				D(2,:,1) = D(2,:,1)/dr
			end if 
		else
			D(2,:,1) = D(2,:,1)/dr
		end if
		!!!!!!!!!!! y(1) !!!!!!!!!!
		D(1,:,2) = cross_product(D(1,:,3),D(1,:,1))
		!!!!!!!!!!! y(2) !!!!!!!!!!
		D(2,:,2) = cross_product(D(2,:,3),D(2,:,1))
		
		!!!!!!!!!!!-----calculating vz-----!!!!!!!!!!!!!!
		h1 = fr%molecule(m)%h1%velocity - fr%molecule(m)%o%velocity
		vz(1)= dot_product(h1,D(1,:,3))
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		h2 = fr%molecule(m)%h2%velocity - fr%molecule(m)%o%velocity
		vz(2)= dot_product(h2,D(2,:,3)) 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		! hydroxyls being evaluated -> only first "water molecule" bond
		if(fr%is_hydroxyl()) then
			b = 1
		else
			b = 2
		end if

		! M_R is calculated
		do n=1,b     ! do on N° oh_bond (=2) X mol
			do i=1,3     ! do on x,y,z (oh_ref)
				MZ(m,time_index)=MZ(m,time_index)+( D(n,bd%R,i)*dMdRz(i)*vz(n) )
			end do
		end do
		
		! A_PQ is calculated
		do n=1,b     ! do on N° oh_bond (=2) X mol
			do i=1,3     ! do on x,y,z (oh_ref)
				do j=1,3     ! do on x,y,z (oh_ref)
					AXX(m,time_index)=AXX(m,time_index)+( D(n,bd%P,i)*dAdRz(i,j)*D(n,bd%Q,j)*vz(n) )
				end do
			end do
		end do
		
		! averages...
		do n = 0,min(t,bd%get_maxlag())
			avMz(m) = avMz(m) + Mz(m,n)/dble(min(t+1,bd%get_maxlag()+1))
			avAxx(m) = avAxx(m) + Axx(m,n)/dble(min(t+1,bd%get_maxlag()+1))
		end do
	  
	end do
	!$OMP END DO
	!$OMP END PARALLEL
end subroutine fill_A_M

subroutine fill_binder(channel, time_index)
	integer*8, intent(in) :: time_index
	integer, intent(IN) :: channel
	integer*1, dimension(bd%NO) :: binder 
	do i=1,bd%NO,2
		read(channel, iostat = ierr) binder(i)
		if(ierr .ne. 0) then
			print*, "ERROR: 'binder.bin' does not have enough data inside"
			stop
		end if
		if((i+1) <= bd%no) then
			binder(i+1) = binder(i) .and. 15 ! lower nibble
		end if
		binder(i) = RSHIFT(binder(i),4) ! upper nibble
	end do
	
	do i = 1, bd%NO
		binder_in_time(i, time_index) = .false.
		do j = 1, size(layers_of_interest)
			if(binder(i) == layers_of_interest(j)) then
				binder_in_time(i, time_index) = .true.
				exit
			end if
		end do
	end do
end subroutine fill_binder

function s_to_HHHMMSS(seconds)
	character(14) :: s_to_HHHMMSS
	real*8, intent(IN) :: seconds
	integer :: hours, minutes, sec
	
	hours = seconds/3600
	minutes = (seconds-hours*3600)/60
	sec = seconds - hours*3600 - minutes*60
	
	write(s_to_HHHMMSS,'(I5.1,A1,I2.2,A1,I2.2)'), hours, ":", minutes, ":", sec
	s_to_HHHMMSS = adjustl(s_to_HHHMMSS)
end function s_to_HHHMMSS

subroutine corr_AMSUM()
! todo AMSUM
	A(nt1) = 0
	MM(nt1) = 0
	do m = 1, bd%NO
		if(binder_in_time(m,nt1)) then
			A(nt1) = A(nt1) + Axx(m,nt1) - AvAxx(m)
			MM(nt1) = MM(nt1) + Mz(m,nt1) - AvMz(m)
		end if
	end do

	!$OMP PARALLEL DEFAULT(SHARED)
	!$OMP DO PRIVATE(timelag, nt0)
	do timelag = 0, min(t,bd%get_maxlag())
		nt0 = mod(t-timelag,corrlen)

		AMSUM_corr(timelag) = AMSUM_corr(timelag) + A(nt1) * MM(nt0)
	end do
	!$OMP END DO
	!$OMP END PARALLEL
end subroutine corr_AMSUM

subroutine corr_self()
	! if skipped self skip frames do the correlation function calculation
	if(mod(t,bd%self_skip) == 0) then 
		!$OMP PARALLEL DEFAULT(SHARED)
		!$OMP DO PRIVATE(timelag, nt0, m)
		! self correlation function
		do timelag = 0, min(t,bd%get_maxlag())
			selfnorm(timelag) = selfnorm(timelag) + 1
			nt0 = mod(t-timelag,corrlen)
			do m = 1, bd%NO
				! if the molecule is not in selected layer continue
				if(.not. (binder_in_time(m,nt0) .or. binder_in_time(m, nt1))) cycle
				
				! molecule M is present in both times nt1 and nt0...
				! full weight of the term
				if(binder_in_time(m,nt0) .and. binder_in_time(m,nt1)) then
					self_corr(timelag) = self_corr(timelag) + (Axx(m, nt1) - avAxx(m)) * (Mz(m, nt0)-avMz(m))
				else
				! molecule M is present only in time nt1 or nt2
				! half weight for the term
					self_corr(timelag) = self_corr(timelag) + 0.5d0 * (Axx(m, nt1) - avAxx(m)) * (Mz(m, nt0)-avMz(m))
				end if
			end do
		end do
		!$OMP END DO
		!$OMP END PARALLEL
	end if
end subroutine corr_self

subroutine corr_cross(neighs)
	logical, intent(in) :: neighs
	integer :: current_neighbor
	! saving oxygen positions...
	do m = 1, bd%NO
		oxygen_position(:,m,nt1) = fr%molecule(m)%O%position
	end do
	! construction of nlist

	! if skipped cross skip frames do the correlation function calculation
	if(mod(t,bd%cross_skip) == 0) then 

		if(neighs .eqv. .true.) then
		
			! CONSTRUCTION OF NLIST
			maxdisplacement = 0 
			secondmaxdisplacement = 0
			maxdisplacementmaxt = 0 
			
			! find maximal displacement over the history memory
			timelag = min(t,bd%get_maxlag())
			nt0 = mod(t-timelag,corrlen)
			!$OMP PARALLEL DEFAULT(SHARED)
			!$OMP DO PRIVATE(m, i, diff, displacement) REDUCTION(max:maxdisplacement, secondmaxdisplacement)
			! this function is written as it is, since I dont trust calling functions inisde parallel region in fortran
			! search for the two greatest displacements to serve as the skin for neighbor list
			do m = 1, bd%NO     
				diff = oxygen_position(:,m,nt0) - oxygen_position(:,m,nt1)
				! pbc check
				do i=1,3
					do while (diff(i) > bd%box_dimensions(i)/2d0)
						diff(i) = diff(i) - bd%box_dimensions(i)
					end do
					
					do while (diff(i) < -bd%box_dimensions(i)/2d0)
						diff(i) = bd%box_dimensions(i) + diff(i)
					end do
				end do

				displacement = norm2(diff)
				
				if(displacement > maxdisplacement) then
					secondmaxdisplacement = maxdisplacement
					maxdisplacement = displacement
				end if
			end do
			!$OMP END DO
			!$OMP END PARALLEL
			
			! the rskin can be set...
			rskin = ccp + maxdisplacement + secondmaxdisplacement
			rskinsq = rskin * rskin

			! establish the NLIST on the level of time t (nt1)
			neighbors = 0
			
			!$OMP PARALLEL DEFAULT(SHARED)
			!$OMP DO PRIVATE(m, n, diff, i, displacement, skip_this)
			do m = 1, bd%NO
				do n = 1, bd%NO
					if(m .eq. n) cycle
					! check if m and n are in rskin distance 
					
					skip_this = .false.                
					! calculate distance between two oxygens in the current frame
					diff = oxygen_position(:,m,nt1) - oxygen_position(:,n,nt1)
					
					! PBC check it is written like this since I don't trust calling subroutines from OMP region...
					do i=1,3
						do while (diff(i) > bd%box_dimensions(i)/2)
							diff(i) = diff(i) - bd%box_dimensions(i)
						end do
					
						do while (diff(i) < -bd%box_dimensions(i)/2)
							diff(i) = bd%box_dimensions(i) + diff(i)
						end do
						
						! if one of the components is bigger than cutoff, this is to potentially save some compute time
						if(diff(i) > rskin) then 
							skip_this = .true.
							exit
						end if
					end do

					! not in cutoff -> go next
					if(skip_this) cycle
					
					! get the square of distance
					displacement = diff(1)**2 + diff(2)**2 + diff(3)**2
					
					! if in cutoff append nlist
					if(displacement <= rskinsq) then
						neighbors(m) = neighbors(m) + 1
						nlist(neighbors(m),m) = n
					end if 
				end do
			end do
			!$OMP END DO
			!$OMP END PARALLEL

			print("(A5,I12,A15,F8.3,A6,f15.11,A6,f15.11,A8,I10,A10,I10)"), "step ", t, " NLIST radius: ", rskin, " dp1: ", maxdisplacement, " dp2: ", secondmaxdisplacement, " pairs: ", sum(neighbors), " avnbcnt: ", sum(neighbors)/bd%no
		end if
	
		! handling of the neighbor list is done

		! calculation of crossterms
		!$OMP PARALLEL DEFAULT(SHARED)
		!$OMP DO PRIVATE(timelag, nt0, m, n, i, diff, displacement, skip_this, auxiliary_variable, current_neighbor)
		do timelag = 0, min(t,bd%get_maxlag())
			crossnorm(timelag) = crossnorm(timelag) + 1
			nt0 = mod(t-timelag,corrlen)
			! go through all molecules
			do m = 1, bd%NO
				! go through all neighbors of molecule m
				do n = 1, neighbors(m)
					! index of nth neighbor of m
					current_neighbor = nlist(n,m)
					! if the molecule m(nt0) and molecule n(nt1) are not in selected layer, go next
					if(.not. (binder_in_time(m,nt0) .or. binder_in_time(current_neighbor, nt1))) cycle
					
					skip_this = .false.                
					! calculate distance between molecule m(nt0) and molecule n(nt1)
					diff = oxygen_position(:,m,nt0) - oxygen_position(:,current_neighbor,nt1)
				
					! PBC check
					do i=1,3
						do while (diff(i) > bd%box_dimensions(i)/2)
							diff(i) = diff(i) - bd%box_dimensions(i)
						end do
					
						do while (diff(i) < -bd%box_dimensions(i)/2)
							diff(i) = bd%box_dimensions(i) + diff(i)
						end do
						
						! if one of the components is bigger than cutoff
						if(diff(i) > ccp) then 
							skip_this = .true.
							exit
						end if
					end do

					! not in cutoff -> go next
					if(skip_this) cycle

					! get the square of distance
					displacement = diff(1)**2 + diff(2)**2 + diff(3)**2
					! if the molecule m(nt0) and molecule n(nt1) are not in cross cut range, go next
					if(displacement > ccpsq) cycle
					
					! the correlation function
					! the weird conditions are the implementation of the switching function

					! if both molecules are in selected layer, add the full term
					if(binder_in_time(m,nt0) .and. binder_in_time(current_neighbor,nt1)) then
						auxiliary_variable = (Axx(current_neighbor, nt1)-avAxx(n)) * (Mz(m, nt0)-avMz(m))
						cross_corr(timelag) = cross_corr(timelag) + auxiliary_variable
						crosscorrM(timelag) = crosscorrM(timelag) + auxiliary_variable
						crosscorrA(timelag) = crosscorrA(timelag) + auxiliary_variable
					else
						! only one molecule of the pair is in the selected layer
						auxiliary_variable = (Axx(current_neighbor, nt1)-avAxx(n)) * (Mz(m, nt0)-avMz(m))
						cross_corr(timelag) = cross_corr(timelag) + 0.5 * auxiliary_variable
						if(binder_in_time(m,nt0)) then
							crosscorrM(timelag) = crosscorrM(timelag) + auxiliary_variable
						else
							crosscorrA(timelag) = crosscorrA(timelag) + auxiliary_variable
						end if
					end if
				end do
			end do
		end do
		!$OMP END DO
		!$OMP END PARALLEL
	end if
end subroutine corr_cross

subroutine corr_hydroxyls()
		! this is quite direct approach to calculate the correlation function, since hydroxyls are allways selected as whole
		! self correlation function
		!$OMP PARALLEL DEFAULT(SHARED)
		!$OMP DO PRIVATE(timelag, nt0, m)
		do timelag = 0, min(t,bd%get_maxlag())
			nt0 = mod(t-timelag,corrlen)
			do m = 1, bd%NO
				hydroxyl_corr(timelag) = hydroxyl_corr(timelag) + (Axx(m, nt1) - avAxx(m)) * (Mz(m, nt0)-avMz(m))
			end do
		end do
		!$OMP END DO
		!$OMP END PARALLEL
end subroutine corr_hydroxyls



end program

	
