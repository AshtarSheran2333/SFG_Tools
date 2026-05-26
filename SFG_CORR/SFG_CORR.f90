program SFG_CORR

!TODO proper modules
use iso_fortran_env
use CORR_UI
use FRAME_READERS
use BOXDATA
use SFG_STRUCTURE
use DXDRZ_DB
use BINDER_FILE
use SFG_CORRELATION

implicit none

#include "utils_error_macros.h"

type(boxdata_type) ::                                       bd

class(frame_reader), allocatable ::                         fr

type(sfg_structure_type) ::                                 struct

type(dXdrz_db_type) ::                                      parameters

type(binder_type) ::                                        binder

type(correlation_function_type) ::                          cf

integer(int64) ::                                           step

integer(int8), dimension(:), allocatable ::                 layer_selection
complex(real64), dimension(:), allocatable ::               spectrum
integer :: i

!evaluate program options
call evaluate_program_options(fr)
allocate(layer_selection(3))
layer_selection(1) = 0
layer_selection(2) = 1
layer_selection(3) = 2

!read BOXDATA
call bd%read_boxdata()

!read structure
error_io_check(struct%read_structure("struct.txt"), "unable to read structure file")

!parameters db
call parameters%init()
error_io_check(parameters%read_from_file("parameters.dat"), "error reading parameters file")

!init binder
error_io_check(binder%init(struct), "unable to init binder")
error_io_check(binder%open_file(), "unable to open binder file")

!CF
call cf%init(struct%groups(1), bd)

!!TODO print recap of the parameters - what groups has been selected... ???
!call recap()

!TODO initialize the correlation function calculation machineries

!TODO timing module would be nice

!TODO the main loop
do step = 1, bd%NSTEP
    !TODO read frame / skip frame
    error_io_check(fr%read_frame(), "unable to read trajectory frame")
    !TODO read binder frame / skip binder frame
    error_io_check(binder%read_frame(), "unable to read binder frame")

    call cf%calculate_step(fr%frame, struct%groups(1), binder%binder_groups(1), parameters, layer_selection, bd) 

    !print progress
    call print_main_loop_progress(step, bd%NSTEP)

    !TODO some way to print backup - sometimes take the correlation data, FFT -> get the convergence series
end do !end of the main loop

call Fourier_transform(cf%correlation_function, bd%DT, bd%DFREQ, bd%FREQ, spectrum, bd%FILTER)
open(84, file = "spectrum.dat")
open(85, file = "corr.dat")
do i = 1, size(spectrum)
    write(84, "(3f16.8)") (i-1)*bd%DFREQ, real(spectrum(i)), imag(spectrum(i))
end do
do i = 1, size(cf%correlation_function)
    write(85, "(2F16.8)") (i-1)*bd%DT, cf%correlation_function(i)
end do
close(84)
close(85)

print*, "DONE"
!TODO finalize

!TODO sort everything below, get rid of it

!contains
!
!subroutine backup()
!    integer :: i, j, f, f_backup, ierr
!    real :: wavenumber
!    complex(kind=8), allocatable, dimension(:,:) :: intermediate_spectrum
!    real(kind=8), allocatable, dimension(:) :: corr_buffer
!    ! 10 times per the trajectory length make spectra from currently calculated correlations
!    ! only the shifted spectra
!    ! only re and im
!    ! for cross skip cross[AM]
!    ! header based on the selected calculation strategies
!    if( mod(t+1, bd%nstep/10) == 0 ) then
!        ! number of frequency points
!        f = bd%FREQ/bd%DFREQ
!        ! allocate buffer for the spectra
!        allocate(intermediate_spectrum(4, 0:f))
!        ! make the buffer zero
!        intermediate_spectrum = 0
!        ! allocate buffer for the normalized corr function
!        allocate(corr_buffer(0:bd%get_maxlag()))
!
!        ! open the backup file
!        if( (t+1)/(bd%nstep/10) == 1 ) then
!            ! open the file in write mode when doing the first backup
!            open(newunit = f_backup, file=trim(output)//"-backup.dat", recl = 160, iostat = ierr)
!            if(ierr .ne. 0) then
!                print*, "error opening", trim(output)//"-backup.dat"
!                return
!            end if
!        else
!            ! open the file in append mode
!            open(newunit = f_backup, file=trim(output)//"-backup.dat", access = 'append', recl = 160, iostat = ierr)
!            if(ierr .ne. 0) then
!                print*, "error opening", trim(output)//"-backup.dat"
!                return
!            end if
!        end if
!        
!        ! calculate FT to the buffers
!        
!        ! write the file
!        
!        write(f_backup,*) "# timestep: ", t+1, " time [ps]: ", (t+1)*bd%dt/1000
!        do i = 0, f
!            wavenumber = dble(i)*bd%DFREQ ! cm^-1
!            write(f_backup,"(f10.2,a)", advance = 'no') wavenumber, " "
!            ! self
!            if(iand(strategy_flags, sf_self) == sf_self) then 
!                if (i .eq. 0) then
!                    ! convert units and normalize the corr function
!                    do j = 0, bd%get_maxlag()
!                        corr_buffer(j) = self_corr(j)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*max(1,selfnorm(j)))
!                    end do
!                    ! calculate the spectrum
!                    call calc_spectrum(corr_buffer, intermediate_spectrum(1,:))
!                    intermediate_spectrum(1,:) = intermediate_spectrum(1,:) - intermediate_spectrum(1,f)
!                end if
!                write(f_backup,"(f16.8,a,f16.8,a)", advance = 'no') real(intermediate_spectrum(1,i)), " ", imag(intermediate_spectrum(1,i)), " "
!            else
!                write(f_backup,"(2I3)", advance = 'no') 0, 0
!            end if
!            ! cross
!            if(iand(strategy_flags, sf_cross) == sf_cross) then 
!                if (i .eq. 0) then
!                    ! convert units and normalize the corr function
!                    do j = 0, bd%get_maxlag()
!                        corr_buffer(j) = cross_corr(j)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*max(1,crossnorm(j)))
!                    end do
!                    ! calculate the spectrum
!                    call calc_spectrum(corr_buffer, intermediate_spectrum(2,:))
!                    intermediate_spectrum(2,:) = intermediate_spectrum(2,:) - intermediate_spectrum(2,f)
!                end if
!                write(f_backup,"(f16.8,a,f16.8,a)", advance = 'no') real(intermediate_spectrum(2,i)), " ", imag(intermediate_spectrum(2,i)), " "
!            else
!                write(f_backup,"(2I3)", advance = 'no') 0, 0
!            end if
!            ! AMSUM
!            if(iand(strategy_flags, sf_AMSUM) == sf_AMSUM) then
!                if (i .eq. 0) then
!                    ! convert units and normalize the corr function
!                    do j = 0, bd%get_maxlag()
!                        corr_buffer(j) = AMSUM_corr(j)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*(t+1-j))
!                    end do
!                    ! calculate the spectrum
!                    call calc_spectrum(corr_buffer, intermediate_spectrum(3,:))
!                    intermediate_spectrum(3,:) = intermediate_spectrum(3,:) - intermediate_spectrum(3,f)
!                end if
!                write(f_backup,"(f16.8,a,f16.8,a)", advance = 'no') real(intermediate_spectrum(3,i)), " ", imag(intermediate_spectrum(3,i)), " "
!            else
!                write(f_backup,"(2I3)", advance = 'no') 0, 0
!            end if
!            ! hydroxyls
!            if(iand(strategy_flags, sf_hydroxyls) == sf_hydroxyls) then
!                if (i .eq. 0) then
!                    ! convert units and normalize the corr function
!                    do j = 0, bd%get_maxlag()
!                        corr_buffer(j) = hydroxyl_corr(j)*(debye_to_ea*electron_to_coulomb)/(bd%box_dimensions(1)*bd%box_dimensions(2)*10d0*(t+1-j))
!                    end do
!                    ! calculate the spectrum
!                    call calc_spectrum(corr_buffer, intermediate_spectrum(4,:))
!                    intermediate_spectrum(4,:) = intermediate_spectrum(4,:) - intermediate_spectrum(4,f)
!                end if
!                write(f_backup,"(f16.8,a,f16.8)", advance = 'no') real(intermediate_spectrum(4,i)), " ", imag(intermediate_spectrum(4,i))
!            else
!                write(f_backup,"(2I3)", advance = 'no') 0, 0
!            end if
!            write(f_backup,*) ""
!        end do
!        ! end the block with extra two empty lines
!        write(f_backup,*) ""
!        write(f_backup,*) ""
!
!        print*, "BACKUP", (t+1)/(bd%nstep/10)
!        
!        ! close the file    
!        close(f_backup)
!        ! deallocate buffer for the spectra
!        deallocate(intermediate_spectrum)
!    end if
!end subroutine backup
!    
!! todo this is ugly...
!subroutine calc_spectrum(in_corr, out_spectrum)
!    real(kind=8), intent(IN), dimension(:) :: in_corr
!    complex(kind=8), intent(INOUT), dimension(:) :: out_spectrum
!    real :: wavenumber, omega, beta, time
!    complex(kind=8) :: eiomegat, eiomegat1
!    integer(kind=8) :: freq, lag, nf
!    
!    out_spectrum = 0
!    
!    nf = bd%FREQ/bd%DFREQ
!    beta = 1/(kb*bd%TEMPERATURE)
!
!    do freq = 0, nf
!        wavenumber = dble(freq)*bd%DFREQ ! cm^-1
!        omega = 2.d0*pi*c*wavenumber ! s^-1
!        
!        ! trapezoidal integration
!        do lag = 1, min(t+1,bd%get_maxlag()+1)
!            time=bd%DT*dble(lag-1)*1.d-15 ! s
!            eiomegat = dcmplx(dcos(dble(omega*time)), dsin(dble(omega*time)))
!            time=bd%DT*dble(lag)*1.d-15 ! s
!            eiomegat1 = dcmplx(dcos(dble(omega*time)), dsin(dble(omega*time)))
!            if(lag < min(t+1,bd%get_maxlag()+1)) then
!                out_spectrum(freq+1) = out_spectrum(freq+1) + 0.5d0 * ( eiomegat * in_corr(lag) * new_filter(lag-1, bd%filter, bd%dt) + eiomegat1 * in_corr(lag+1) * new_filter(lag, bd%filter, bd%dt) ) * bd%DT
!            else ! add one artificial point to be 0 for the trap rule...
!                out_spectrum(freq+1) = out_spectrum(freq+1) + 0.5d0 * eiomegat * in_corr(lag) * new_filter(lag-1, bd%filter, bd%dt) * bd%DT
!            end if
!        end do
!        
!        out_spectrum(freq+1) = dcmplx(0,-1) * beta * out_spectrum(freq+1) / omega
!    end do
!    
!
!end subroutine calc_spectrum
!    
!! make_spectrum(arg1, arg2)
!! arg1 - correlation function array
!! arg2 - prefix of the spectrum file
!! makes Fourier transform of correlation function, multiplies it by 1/kbT, and saves the spectrum to the file
!! todo this is ugly
!subroutine make_spectrum(corr, name)
!    real(8), dimension(:), intent(IN) ::                    corr
!    
!    character(len=*), intent(IN) ::                         name
!    
!    complex(8), allocatable, dimension(:) ::                spectrum
!    
!    real(8) ::                                              wavenumber
!    
!    integer(8) ::                                           f,&
!                                                            f_spectrum,&
!                                                            f_shiftedspectrum
!    
!    f = bd%FREQ/bd%DFREQ
!
!    allocate(spectrum(f+1))
!    spectrum = 0
!    
!    call calc_spectrum(corr, spectrum)
!    
!    print*,""
!    print*, name//"-spectrum.dat"
!    open(newunit = f_spectrum, file=name//"-spectrum.dat", recl=120, iostat = ierr)
!    if(ierr .ne. 0) then
!        print"(a,a,a)", "ERROR: unable to create file '", name//"-spectrum.dat", "'"
!        stop
!    end if
!    write(f_spectrum,"(A,F8.3)"), "#FILTER PARAMETER: ", bd%filter ! in ps
!    write(f_spectrum,"(A)") "#wavelength cm-1, re, im, abs, phase deg"
!    print*, name//"-shiftedspectrum.dat"
!    print*, name//"-shiftedspectrum.dat is shifted by:"
!    print*, "(Re, Im)", -spectrum(f)
!    print*, "which is value of the spectrum at frequency ", f, "cm^-1"
!    open(newunit = f_shiftedspectrum, file=name//"-shiftedspectrum.dat", recl=120, iostat = ierr)
!    if(ierr .ne. 0) then
!        print"(a,a,a)", "ERROR: unable to create file '", name//"-shiftedspectrum.dat", "'"
!        stop
!    end if
!    write(f_shiftedspectrum,"(A,F8.3)"), "#FILTER PARAMETER: ", bd%filter ! in ps
!    write(f_shiftedspectrum,"(A)") "#wavelength cm-1, re, im, abs, phase deg"
!    
!    do i=1,f+1
!        ! according to Khatib equation (3) gets the second order susceptibility
!        wavenumber=dble(i-1)*bd%DFREQ       ! cm-1
!    
!        write(f_spectrum,*) int(wavenumber),&
!                    real(spectrum(i)),&
!                    imag(spectrum(i)),&
!                    abs(spectrum(i)),&
!                    datan(imag(spectrum(i))/real(spectrum(i)))*180/pi
!     
!        ! the spectrum is shifted by the value at frequency $FREQ, sine there should be no vibrations
!        write(f_shiftedspectrum,*) wavenumber,&
!                    real(spectrum(i)-spectrum(f+1)),&
!                    imag(spectrum(i)-spectrum(f+1)),&
!                    abs(spectrum(i)-spectrum(f+1)),&
!                    datan(imag(spectrum(i)-spectrum(f+1))/real(spectrum(i)-spectrum(f+1)))*180/pi
!    enddo
!    
!    close(f_spectrum)
!    close(f_shiftedspectrum)
!    
!    deallocate(spectrum)
!end subroutine make_spectrum
!
!real*8 function new_filter(int_time, real_filter, dt)
!    integer(kind=8), intent(IN) :: int_time ! steps
!    real(kind=8), intent(IN) :: dt, real_filter !fs, ps
!
!    new_filter = dexp(-(int_time*dt/(1000*real_filter))**2)
!
!end function new_filter
!
!real*8 function filter(time,tau)
!    integer*8, intent(IN) ::                                time
!    real*8, intent(IN) ::                                   tau
!    
!    filter = dexp(-(dble(time)/dble(tau))**2)
!end
!   
!logical function in_range(molA, molB, timeA, timeB, cross_corr_parameter, cross_corr_parameter2)
!    integer*8, intent(IN) :: molA, molB, timeA, timeB
!    real*8, intent(IN) :: cross_corr_parameter, cross_corr_parameter2
!    real*8, dimension(3) :: posDiff
!    real*8 :: distance
!    integer*1 :: i
!    
!    in_range = .FALSE.
!    ! calculate distance between two oxygens
!    posDiff = oxygen_position(:,molA,timeA) - oxygen_position(:,molB,timeB)
!    
!    ! PBC check
!    do i=1,3
!        if (posDiff(i) > bd%box_dimensions(i)/2) then
!            posDiff(i) = posDiff(i) - bd%box_dimensions(i)
!        else if (posDiff(i) < -bd%box_dimensions(i)/2) then
!            posDiff(i) = bd%box_dimensions(i) + posDiff(i)
!        end if
!        if(posDiff(i) > cross_corr_parameter) return
!    end do
!    
!    distance = posDiff(1)**2 + posDiff(2)**2 + posDiff(3)**2
!    
!    if(distance <= cross_corr_parameter2) in_range = .TRUE.
!end function in_range
!
!subroutine recap()
!    print*, ""
!    print*, "Recap of calculation settings:"
!    print*, ""
!    if(output .ne. "") then
!        print*, "Output file prefix: ", trim(output)
!        print*, ""
!    else
!        print*, "Output file prefix: ", "-no prefix selected-"
!        print*, ""
!    end if
!    
!    if(iand(strategy_flags, sf_self) .eq. sf_self) then ! self
!        print*, "self correlation terms will be calculated"
!        print*, "self correlation function"
!        print*, "will be stored in file: <output file prefix>-selfterms.dat"
!        print*, "self correlation spectrum"
!        print*, "will be stored in files: <output file prefix>-self-*.dat"
!        print*, ""
!    end if
!
!    if(iand(strategy_flags, sf_cross) .eq. sf_cross) then ! cross
!        print*, "cross correlation terms will be calculated"
!        if(enable_nlist .eqv. .true.) then
!            print*, "neighborlist will be used"
!        else
!            print*, "neighbor list will not be used"
!        end if
!        print*, "r_cut: ", ccp, "Angstrom"
!        print*, "cross correlation function"
!        print*, "will be stored in file: <output file prefix>-crossterms.dat"
!        print*, "cross correlation spectrum"
!        print*, "will be stored in files: <output file prefix>-cross-*.dat"
!        print*, ""
!    end if
!
!    if(iand(strategy_flags, sf_AMSUM) .eq. sf_AMSUM) then ! AMSUM
!        print*, "corelation of whole selected layer summed properties will be calculated"
!        print*, "the correlation function"
!        print*, "will be stored in file: <output file prefix>-AMSUMterms.dat"
!        print*, "the spectrum"
!        print*, "will be stored in file: <output file prefix>-AMSUMterms.dat"
!        print*, ""
!    end if
!
!    print*, "dt is set to [fs]: ", bd%DT
!    print*, "correlation function length [ps]:", bd%get_maxlag()*bd%dt/1000
!    print*, "temperature is set to [K]: ", bd%TEMPERATURE
!    print*, "Filter parameter is set to [ps]: ", bd%filter ! ps
!    print*, "" ! dM/dr dA/dr
!    print*, "dM_x/dr_z: ", dMdRz(1)
!    print*, "dM_x/dr_z: ", dMdRz(2)
!    print*, "dM_x/dr_z: ", dMdRz(3)
!    print*, ""
!    print*, "dA_xx/dr_z: ", dAdRz(1,1)
!    print*, "dA_yy/dr_z: ", dAdRz(2,2)
!    print*, "dA_zz/dr_z: ", dAdRz(3,3)
!    print*, "dA_xy/dr_z: ", dAdRz(1,2)
!    print*, "dA_xz/dr_z: ", dAdRz(1,3)
!    print*, "dA_yz/dr_z: ", dAdRz(2,3)
!    
!    print*, "_______________________________________________________________________________"
!end subroutine recap
!
!subroutine allocate_memory()
!    ! allocate and set general variables
!    allocate(Axx(bd%NO, 0:bd%get_maxlag()))
!    allocate(Mz(bd%NO, 0:bd%get_maxlag()))
!    allocate(avMz(bd%NO))
!    allocate(avAXX(bd%NO))
!    axx = 0
!    mz = 0
!    avMz=0
!    avAXX=0
!
!    ! strategy hydroxyls
!    if(strategy_flags < 0) then 
!        allocate(hydroxyl_corr(0:bd%get_maxlag()))
!        hydroxyl_corr = 0
!    else
!        ! if not dealing with hydroxyls, binder is needed
!        allocate(binder_in_time(BD%NO,0:bd%get_maxlag()))
!
!        ! strategy self
!        if(iand(strategy_flags, sf_self) == sf_self) then
!            allocate(self_corr(0:bd%get_maxlag()))
!            allocate(selfnorm(0:bd%get_maxlag()))
!            selfnorm = 0
!            self_corr = 0
!        end if
!        
!        ! strategy cross
!        if(iand(strategy_flags, sf_cross) == sf_cross) then
!            allocate(cross_corr(0:bd%get_maxlag()))
!            allocate(crosscorrM(0:bd%get_maxlag()))
!            allocate(crosscorrA(0:bd%get_maxlag()))
!            crosscorrM = 0
!            crosscorrA = 0
!            allocate(oxygen_position(3,bd%NO,0:bd%get_maxlag()))
!            cross_corr = 0
!            ! NLIST
!            allocate(nlist(bd%NO,bd%NO))
!            nlist = 0
!            ! CAREFUL first dimension can be shrinked...
!            ! this did not cause any problems since our systems are < 10k water molecules
!            ! this could be done easily in c++ utilizing std::vector...
!            ! I dont want to implement dynamic array RN since this is just experimental feature
!            ! I can theoretically add another calculation strategy for cross_nlist... 
!            ! to save all this space... and some computation costs
!            allocate(neighbors(bd%NO))
!            neighbors = 0
!            allocate(crossnorm(0:bd%get_maxlag()))
!            crossnorm = 0
!            
!            if(.not. enable_nlist) then
!                ! just fill the neighbor list with all possible combinantions
!                do m = 1, bd%NO
!                    do n = 1, bd%NO
!                        if(m .eq. n) cycle
!                        neighbors(m) = neighbors(m) + 1
!                        nlist(neighbors(m),m) = n
!                    end do
!                end do
!            end if
!        end if
!        
!        ! strategy AMSUM
!        if(iand(strategy_flags, sf_AMSUM) == sf_AMSUM) then
!            allocate(A(0:bd%get_maxlag()))
!            A = 0
!            allocate(MM(0:bd%get_maxlag()))
!            MM = 0
!            allocate(AMSUM_corr(0:bd%get_maxlag()))
!            AMSUM_corr = 0
!        end if
!    end if
!    
!end subroutine allocate_memory
!
!subroutine fill_A_M(time_index)
!    integer*8, intent(IN) :: time_index
!    integer :: m,n,i,j,b=2
!    real*8 :: dr, scalar
!    real*8, dimension(3) :: h1, h2
!    real*8, dimension(2,3,3) :: D
!    real*8, dimension(2) :: vz
!
!    Axx(:,time_index) = 0
!    Mz(:,time_index) = 0
!    
!    avMz=0.d0
!    avAXX=0.d0
!
!    ! print*, "evaluating D matrices and v(z) for each mol"
!    !$OMP PARALLEL DEFAULT(SHARED)
!    !$OMP DO PRIVATE(m,n,i,j,h1,h2,dr,scalar,D,vz)
!    do m=1,bd%NO
!        
!        !!!!!!!!!!! z(1) !!!!!!!!!!
!        h1 = fr%molecule(m)%h1%position - fr%molecule(m)%o%position
!        ! pbc check
!        do i=1,3
!            do while (h1(i) > bd%box_dimensions(i)/2d0)
!                h1(i) = h1(i) - bd%box_dimensions(i)
!            end do
!            
!            do while (h1(i) < -bd%box_dimensions(i)/2d0)
!                h1(i) = bd%box_dimensions(i) + h1(i)
!            end do
!        end do
!        dr = norm2(h1)
!        D(1,:,3) = h1/dr
!        !!!!!!!!!!! z(2) !!!!!!!!!!
!        h2 = fr%molecule(m)%h2%position - fr%molecule(m)%o%position
!        ! pbc check
!        do i=1,3
!            do while (h2(i) > bd%box_dimensions(i)/2d0)
!                h2(i) = h2(i) - bd%box_dimensions(i)
!            end do
!            
!            do while (h2(i) < -bd%box_dimensions(i)/2d0)
!                h2(i) = bd%box_dimensions(i) + h2(i)
!            end do
!        end do
!        dr = norm2(h2)
!        D(2,:,3) = h2/dr
!        !!!!!!!!!!! x(1) !!!!!!!!!!
!        scalar = dot_product(h2,D(1,:,3))
!        D(1,:,1) = scalar*D(1,:,3)-h2
!        dr = norm2(D(1,:,1))
!        if(dr == 0) then 
!            if((D(1,1,3) == 0) .and. (D(1,2,3) == 0)) then
!                ! be careful, in this case we would divide by zero!
!                ! X axis if lab Z == Z else -X axis
!                D(1,:,1) = (/D(1,3,3),0d0,0d0/)
!            else
!                ! if hoh angle is 180deg (should be impossible in case of water molecules - but useful in case of surface hydroxyls)
!                ! chose X perpendicular to Z
!                D(1,:,1) = (/D(1,2,3),-D(1,1,3),0d0/)
!                dr = norm2(D(1,:,1))
!                D(1,:,1) = D(1,:,1)/dr
!            end if
!        else
!            D(1,:,1) = D(1,:,1)/dr
!        end if
!        !!!!!!!!!!! x(2) !!!!!!!!!!
!        scalar = dot_product(h1,D(2,:,3))
!        D(2,:,1) = scalar*D(2,:,3)-h1
!        dr = norm2(D(2,:,1))
!        if(dr == 0) then 
!            if((D(2,1,3) == 0) .and. (D(2,2,3) == 0)) then
!                ! be careful, in this case we would divide by zero!
!                ! X axis if lab Z == Z else -X axis
!                D(2,:,1) = (/D(2,3,3),0d0,0d0/)
!            else
!                ! if hoh angle is 180deg (should be impossible in case of water molecules - but useful in case of surface hydroxyls)
!                ! chose X perpendicular to Z
!                D(2,:,1) = (/D(2,2,3),-D(2,1,3),0d0/)
!                dr = norm2(D(2,:,1))
!                D(2,:,1) = D(2,:,1)/dr
!            end if 
!        else
!            D(2,:,1) = D(2,:,1)/dr
!        end if
!        !!!!!!!!!!! y(1) !!!!!!!!!!
!        D(1,:,2) = cross_product(D(1,:,3),D(1,:,1))
!        !!!!!!!!!!! y(2) !!!!!!!!!!
!        D(2,:,2) = cross_product(D(2,:,3),D(2,:,1))
!        
!        !!!!!!!!!!!-----calculating vz-----!!!!!!!!!!!!!!
!        h1 = fr%molecule(m)%h1%velocity - fr%molecule(m)%o%velocity
!        vz(1)= dot_product(h1,D(1,:,3))
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        h2 = fr%molecule(m)%h2%velocity - fr%molecule(m)%o%velocity
!        vz(2)= dot_product(h2,D(2,:,3)) 
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        
!        ! hydroxyls being evaluated -> only first "water molecule" bond
!        if(fr%is_hydroxyl()) then
!            b = 1
!        else
!            b = 2
!        end if
!
!        ! M_R is calculated
!        do n=1,b     ! do on N° oh_bond (=2) X mol
!            do i=1,3     ! do on x,y,z (oh_ref)
!                MZ(m,time_index)=MZ(m,time_index)+( D(n,bd%R,i)*dMdRz(i)*vz(n) )
!            end do
!        end do
!        
!        ! A_PQ is calculated
!        do n=1,b     ! do on N° oh_bond (=2) X mol
!            do i=1,3     ! do on x,y,z (oh_ref)
!                do j=1,3     ! do on x,y,z (oh_ref)
!                    AXX(m,time_index)=AXX(m,time_index)+( D(n,bd%P,i)*dAdRz(i,j)*D(n,bd%Q,j)*vz(n) )
!                end do
!            end do
!        end do
!        
!        ! averages...
!        do n = 0,min(t,bd%get_maxlag())
!            avMz(m) = avMz(m) + Mz(m,n)/dble(min(t+1,bd%get_maxlag()+1))
!            avAxx(m) = avAxx(m) + Axx(m,n)/dble(min(t+1,bd%get_maxlag()+1))
!        end do
!      
!    end do
!    !$OMP END DO
!    !$OMP END PARALLEL
!end subroutine fill_A_M
!
!subroutine corr_AMSUM()
!! todo AMSUM
!    A(nt1) = 0
!    MM(nt1) = 0
!    do m = 1, bd%NO
!        if(binder_in_time(m,nt1)) then
!            A(nt1) = A(nt1) + Axx(m,nt1) - AvAxx(m)
!            MM(nt1) = MM(nt1) + Mz(m,nt1) - AvMz(m)
!        end if
!    end do
!
!    !$OMP PARALLEL DEFAULT(SHARED)
!    !$OMP DO PRIVATE(timelag, nt0)
!    do timelag = 0, min(t,bd%get_maxlag())
!        nt0 = mod(t-timelag,corrlen)
!
!        AMSUM_corr(timelag) = AMSUM_corr(timelag) + A(nt1) * MM(nt0)
!    end do
!    !$OMP END DO
!    !$OMP END PARALLEL
!end subroutine corr_AMSUM
!
!subroutine corr_self()
!    ! if skipped self skip frames do the correlation function calculation
!    if(mod(t,bd%self_skip) == 0) then 
!        !$OMP PARALLEL DEFAULT(SHARED)
!        !$OMP DO PRIVATE(timelag, nt0, m)
!        ! self correlation function
!        do timelag = 0, min(t,bd%get_maxlag())
!            selfnorm(timelag) = selfnorm(timelag) + 1
!            nt0 = mod(t-timelag,corrlen)
!            do m = 1, bd%NO
!                ! if the molecule is not in selected layer continue
!                if(.not. (binder_in_time(m,nt0) .or. binder_in_time(m, nt1))) cycle
!                
!                ! molecule M is present in both times nt1 and nt0...
!                ! full weight of the term
!                if(binder_in_time(m,nt0) .and. binder_in_time(m,nt1)) then
!                    self_corr(timelag) = self_corr(timelag) + (Axx(m, nt1) - avAxx(m)) * (Mz(m, nt0)-avMz(m))
!                else
!                ! molecule M is present only in time nt1 or nt2
!                ! half weight for the term
!                    self_corr(timelag) = self_corr(timelag) + 0.5d0 * (Axx(m, nt1) - avAxx(m)) * (Mz(m, nt0)-avMz(m))
!                end if
!            end do
!        end do
!        !$OMP END DO
!        !$OMP END PARALLEL
!    end if
!end subroutine corr_self
!
!subroutine corr_cross(neighs)
!    logical, intent(in) :: neighs
!    integer :: current_neighbor
!    ! saving oxygen positions...
!    do m = 1, bd%NO
!        oxygen_position(:,m,nt1) = fr%molecule(m)%O%position
!    end do
!    ! construction of nlist
!
!    ! if skipped cross skip frames do the correlation function calculation
!    if(mod(t,bd%cross_skip) == 0) then 
!
!        if(neighs .eqv. .true.) then
!        
!            ! CONSTRUCTION OF NLIST
!            maxdisplacement = 0 
!            secondmaxdisplacement = 0
!            maxdisplacementmaxt = 0 
!            
!            ! find maximal displacement over the history memory
!            timelag = min(t,bd%get_maxlag())
!            nt0 = mod(t-timelag,corrlen)
!            !$OMP PARALLEL DEFAULT(SHARED)
!            !$OMP DO PRIVATE(m, i, diff, displacement) REDUCTION(max:maxdisplacement, secondmaxdisplacement)
!            ! this function is written as it is, since I dont trust calling functions inisde parallel region in fortran
!            ! search for the two greatest displacements to serve as the skin for neighbor list
!            do m = 1, bd%NO     
!                diff = oxygen_position(:,m,nt0) - oxygen_position(:,m,nt1)
!                ! pbc check
!                do i=1,3
!                    do while (diff(i) > bd%box_dimensions(i)/2d0)
!                        diff(i) = diff(i) - bd%box_dimensions(i)
!                    end do
!                    
!                    do while (diff(i) < -bd%box_dimensions(i)/2d0)
!                        diff(i) = bd%box_dimensions(i) + diff(i)
!                    end do
!                end do
!
!                displacement = norm2(diff)
!                
!                if(displacement > maxdisplacement) then
!                    secondmaxdisplacement = maxdisplacement
!                    maxdisplacement = displacement
!                end if
!            end do
!            !$OMP END DO
!            !$OMP END PARALLEL
!            
!            ! the rskin can be set...
!            rskin = ccp + maxdisplacement + secondmaxdisplacement
!            rskinsq = rskin * rskin
!
!            ! establish the NLIST on the level of time t (nt1)
!            neighbors = 0
!            
!            !$OMP PARALLEL DEFAULT(SHARED)
!            !$OMP DO PRIVATE(m, n, diff, i, displacement, skip_this)
!            do m = 1, bd%NO
!                do n = 1, bd%NO
!                    if(m .eq. n) cycle
!                    ! check if m and n are in rskin distance 
!                    
!                    skip_this = .false.                
!                    ! calculate distance between two oxygens in the current frame
!                    diff = oxygen_position(:,m,nt1) - oxygen_position(:,n,nt1)
!                    
!                    ! PBC check it is written like this since I don't trust calling subroutines from OMP region...
!                    do i=1,3
!                        do while (diff(i) > bd%box_dimensions(i)/2)
!                            diff(i) = diff(i) - bd%box_dimensions(i)
!                        end do
!                    
!                        do while (diff(i) < -bd%box_dimensions(i)/2)
!                            diff(i) = bd%box_dimensions(i) + diff(i)
!                        end do
!                        
!                        ! if one of the components is bigger than cutoff, this is to potentially save some compute time
!                        if(diff(i) > rskin) then 
!                            skip_this = .true.
!                            exit
!                        end if
!                    end do
!
!                    ! not in cutoff -> go next
!                    if(skip_this) cycle
!                    
!                    ! get the square of distance
!                    displacement = diff(1)**2 + diff(2)**2 + diff(3)**2
!                    
!                    ! if in cutoff append nlist
!                    if(displacement <= rskinsq) then
!                        neighbors(m) = neighbors(m) + 1
!                        nlist(neighbors(m),m) = n
!                    end if 
!                end do
!            end do
!            !$OMP END DO
!            !$OMP END PARALLEL
!
!            print("(A5,I12,A15,F8.3,A6,f15.11,A6,f15.11,A8,I10,A10,I10)"), "step ", t, " NLIST radius: ", rskin, " dp1: ", maxdisplacement, " dp2: ", secondmaxdisplacement, " pairs: ", sum(neighbors), " avnbcnt: ", sum(neighbors)/bd%no
!        end if
!    
!        ! handling of the neighbor list is done
!
!        ! calculation of crossterms
!        !$OMP PARALLEL DEFAULT(SHARED)
!        !$OMP DO PRIVATE(timelag, nt0, m, n, i, diff, displacement, skip_this, auxiliary_variable, current_neighbor)
!        do timelag = 0, min(t,bd%get_maxlag())
!            crossnorm(timelag) = crossnorm(timelag) + 1
!            nt0 = mod(t-timelag,corrlen)
!            ! go through all molecules
!            do m = 1, bd%NO
!                ! go through all neighbors of molecule m
!                do n = 1, neighbors(m)
!                    ! index of nth neighbor of m
!                    current_neighbor = nlist(n,m)
!                    ! if the molecule m(nt0) and molecule n(nt1) are not in selected layer, go next
!                    if(.not. (binder_in_time(m,nt0) .or. binder_in_time(current_neighbor, nt1))) cycle
!                    
!                    skip_this = .false.                
!                    ! calculate distance between molecule m(nt0) and molecule n(nt1)
!                    diff = oxygen_position(:,m,nt0) - oxygen_position(:,current_neighbor,nt1)
!                
!                    ! PBC check
!                    do i=1,3
!                        do while (diff(i) > bd%box_dimensions(i)/2)
!                            diff(i) = diff(i) - bd%box_dimensions(i)
!                        end do
!                    
!                        do while (diff(i) < -bd%box_dimensions(i)/2)
!                            diff(i) = bd%box_dimensions(i) + diff(i)
!                        end do
!                        
!                        ! if one of the components is bigger than cutoff
!                        if(diff(i) > ccp) then 
!                            skip_this = .true.
!                            exit
!                        end if
!                    end do
!
!                    ! not in cutoff -> go next
!                    if(skip_this) cycle
!
!                    ! get the square of distance
!                    displacement = diff(1)**2 + diff(2)**2 + diff(3)**2
!                    ! if the molecule m(nt0) and molecule n(nt1) are not in cross cut range, go next
!                    if(displacement > ccpsq) cycle
!                    
!                    ! the correlation function
!                    ! the weird conditions are the implementation of the switching function
!
!                    ! if both molecules are in selected layer, add the full term
!                    if(binder_in_time(m,nt0) .and. binder_in_time(current_neighbor,nt1)) then
!                        auxiliary_variable = (Axx(current_neighbor, nt1)-avAxx(n)) * (Mz(m, nt0)-avMz(m))
!                        cross_corr(timelag) = cross_corr(timelag) + auxiliary_variable
!                        crosscorrM(timelag) = crosscorrM(timelag) + auxiliary_variable
!                        crosscorrA(timelag) = crosscorrA(timelag) + auxiliary_variable
!                    else
!                        ! only one molecule of the pair is in the selected layer
!                        auxiliary_variable = (Axx(current_neighbor, nt1)-avAxx(n)) * (Mz(m, nt0)-avMz(m))
!                        cross_corr(timelag) = cross_corr(timelag) + 0.5 * auxiliary_variable
!                        if(binder_in_time(m,nt0)) then
!                            crosscorrM(timelag) = crosscorrM(timelag) + auxiliary_variable
!                        else
!                            crosscorrA(timelag) = crosscorrA(timelag) + auxiliary_variable
!                        end if
!                    end if
!                end do
!            end do
!        end do
!        !$OMP END DO
!        !$OMP END PARALLEL
!    end if
!end subroutine corr_cross
!
!subroutine corr_hydroxyls()
!        ! this is quite direct approach to calculate the correlation function, since hydroxyls are allways selected as whole
!        ! self correlation function
!        !$OMP PARALLEL DEFAULT(SHARED)
!        !$OMP DO PRIVATE(timelag, nt0, m)
!        do timelag = 0, min(t,bd%get_maxlag())
!            nt0 = mod(t-timelag,corrlen)
!            do m = 1, bd%NO
!                hydroxyl_corr(timelag) = hydroxyl_corr(timelag) + (Axx(m, nt1) - avAxx(m)) * (Mz(m, nt0)-avMz(m))
!            end do
!        end do
!        !$OMP END DO
!        !$OMP END PARALLEL
!end subroutine corr_hydroxyls

contains


end program

    
