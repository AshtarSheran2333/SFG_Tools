program interf

use SFG_TYPES
use bdata
use SFG_FRAMES
use switches

implicit none

type(boxdata) ::                                            bd

class(frame_reader), pointer ::                             fr
type(trr_frame_reader), allocatable, target ::              trr
type(gro_frame_reader), allocatable, target ::              gro
type(xyz_frame_reader), allocatable, target ::              xyz

real*8, dimension (:,:,:), allocatable ::                   points
!points replaces lx,ly,lz --- coordinates of both interfaces

real*8, dimension(:,:), allocatable ::                      density_function
!holds up all the density functions (distance,density,updensity,downdensity)
real*8, dimension(:), allocatable ::						density_axis

real*8, dimension(3) ::                                     prev_rod_pos,&
															rod_pos,&
															rod_start,&
															diff
!rod_pos replaces l1, l2, l3 --- coordinates of "poking rod"
!diff replaces xdiff, ydiff, zdiff --- stores temporary distance of oxygen atoms from "poking rod"
!rod_start replaces fx, fy, fz --- coordinates of rod starting position

integer*8, dimension (:,:,:), allocatable ::                prev_index

integer*8, dimension(:,:), allocatable ::                   ij_array

integer*8, dimension(2,2) ::                                ij_est

integer*8, dimension(3) ::                                  n_points
!n_points replaces nlx, nly, nlz --- number of divisions in each direction given by dividing boxdimensions by dl

integer, dimension(2) ::                                    uji,&
															dji

integer ::                                                  ierr,&
															f_interface,&
															f_grid_interface,&
															f_interfacebin

real*8 ::                                                   r,&
															p,&
															pdiff,&
															u,&
															urmin,&
															drmin,&
															uzdiff,&
															dzdiff,&
															zcenter,&
															search_radius,&
															prev_pdiff
!r - distance of given oxygen atom from the end of "poking rod"
!p - sum[oxygen atoms] if condition is met - exp(-r**2/(2*E**2))/((2*pi*E**2)**1.5)
!u
!urmin - min of !up! r = ||(xdiff,ydiff,zdiff)||
!drmin - min of !down! r = ||(xdiff,ydiff,zdiff)||
!uzdiff 
!dzdiff
!zcenter - z coordinate approximately in the middle of water slab

integer*8 ::                                                step,&
															d,&
															i,&
															j,&
															k,&
															m,&
															s,&
															nop,&
															temp
!step - iterator through steps
!d - sumcheck of points each frame -> should be allways n_points(1)*n_points(2)*2
!i,j,k -general iterators
!m -iterates through NO
!s -temporary storage for writing step number to the output file

logical ::                                                  found_uper_interface,&
															found_bottom_interface,&
															vmdout = .FALSE.,&
															vmdpbc = .TRUE.,&
															bin_header_done = .FALSE.,&
															nocheck = .false.
!found_uper_interface/min holds up info if given point met conditions to become rod_max/min  

character(len = 256) ::                                     errmsg

character(len=60), allocatable ::                           arg,&
															arg1

real*8, parameter ::                                        pi=dacos(-1.d0),&
															E=2.4,&
															tollerance=0.004,&
															waterDensity=0.03336 !N/A^3

!##########################################EVALUATE PROGRAM SWITCHES###########################################

call evaluate_switches()

!###############################################################################################################

call bd%read_boxdata()

if(.NOT. nocheck) then
	print*, "Checking the contents of trajectory file(s)"
	print*, ""
	do i = 1, bd%NSTEP
		call fr%skip_frame()
	end do
	call fr%rewind_file()
	print*, "contents of trajectory files - OK"
	print*, ""
end if

print*, "Interface parameters recap:"
print*, ""
print*, "Instantaneous surface parameters:"
print*, ""
print*, "resolution of instantaneous surface (x,y,z) is:"
print*, bd%interface_volume_element
print*, ""
print*, "number of skipped frames:", bd%INTERFACE_SKIP
print*, ""
print*, "Density function parameters:"
print*, ""
print*, "Distance from instantaneous surface analyzed (Angstrom):"
print*, bd%density_min_r, "to", bd%density_max_r
print*, ""
print*, "number of bins in 1 A:", bd%density_bin
print*, ""
print*, "width of bin (Angstrom): ", 2*bd%density_tol
print*, ""
print*, "_______________________________________________________________________________"

!##############################################################################################################

if(vmdout) then 
	open (newunit = f_interface, FILE="interface.xyz", iostat = ierr)
	if(ierr .ne. 0) then
		print*, "ERROR: -vmdout is unable to create a file called 'interface.xyz'" 
		stop
	end if
end if

open (newunit = f_grid_interface, FILE="grid_interface.dat", iostat = ierr)
if(ierr .ne. 0) then
	print*, "ERROR: unable to create file called 'grid_interface.dat'" 
	stop
end if
	
if(bd%cancel_interface_calculation == .false.) then
	!if calculation is enabled, create file
	open (newunit = f_interfacebin, FILE = "interface.bin", form = "unformatted", access="stream", convert = 'big_endian', iostat = ierr)
	if(ierr .ne. 0) then
		print*, "ERROR: unable to create file called 'interface.bin'" 
		stop
	end if
else
	!else open
	open (newunit = f_interfacebin, FILE = "interface.bin", form = "unformatted", access="stream", convert = 'big_endian', status = 'old', action = 'read', iostat = ierr)
	if(ierr .ne. 0) then
		print*, "ERROR: unable to open 'interface.bin'" 
		stop
	end if
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! density file preparation
nop=(bd%density_max_r-bd%density_min_r)*bd%density_bin+1

allocate (density_function(nop,3))
allocate (density_axis(nop))
density_function = 0

do i=1,nop
	density_axis(i)=(real(i,8)-1)/bd%density_bin+bd%density_min_r
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!   GRID DEFINITION     !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

n_points = NINT(bd%box_dimensions/bd%interface_volume_element)

bd%interface_volume_element = bd%box_dimensions/n_points

rod_start = -bd%box_dimensions(:)/2.0 + bd%interface_volume_element(:)/2.0

n_points(:) = n_points(:)-1

write(f_grid_interface,*) "========================="
write(f_grid_interface,*) " Npoint | first | last"
do i = 1, 3
	write(f_grid_interface,*) n_points(i), rod_start(i), n_points(i)*bd%interface_volume_element(i)+rod_start(i)
end do
write(f_grid_interface,*) "========================="
close(f_grid_interface)


allocate(points(0:n_points(1),0:n_points(2),2))

allocate(prev_index(0:n_points(1),0:n_points(2),2))
!fills zeros to previous bottom interface indexes
prev_index(:,:,2) = 0
!fills max value to previous top interface index
prev_index(:,:,1) = 10000000
!the values will be rewriten when the frame 1 interface is found and next interface is calculated according to that

do step = 0, bd%NSTEP-1, bd%INTERFACE_SKIP
	
	print "(a,a,i,$)", char(13), "STEP: ", step+1
	
	call fr%read_frame() !reading frame
	
	!estimating center of the water slab
	if (step == 0) then
		print*, ""
		print*, "estimating center of the water slab..."
		zcenter = 0
		do i = 1, bd%NO
			zcenter = zcenter + fr%molecule(i)%o%position(3)
		end do

		if(bd%NO .ne. size(fr%molecule(:))) then
			print*, "number of oxygens in BOXDATA file ", bd%NO, " does not match the number of oxygens in the frame ", size(fr%molecule(:))
			exit
		endif

		zcenter = zcenter/bd%NO !searching center of water slab as center of mass of all the oxygens
		print*, "estimated center of the water slab is at Z = ", zcenter
		print*, "_______________________________________________________________________________"
	end if
	
	if(bd%cancel_interface_calculation .eq. .true.) then
		call read_interface_frame()
	else
		!get the interface
		!poking rod through all grid points
		!$OMP PARALLEL DEFAULT(NONE) SHARED(n_points, bd, prev_index, points, fr, step, rod_start) PRIVATE( i, j, k, m, s, found_uper_interface, found_bottom_interface, rod_pos, pdiff, p, diff, r, prev_pdiff, prev_rod_pos, errmsg )
		!$OMP DO
		do i=0,n_points(1)
	
			do j=0,n_points(2)
			
				found_uper_interface= .false.
				found_bottom_interface= .false.
				
				rod_pos = rod_start
				pdiff = bd%interface_density*waterDensity
	
				!poking rod from bottom
				do k=max(0,prev_index(i,j,2)-bd%INTERFACE_PUSHBACK) ,n_points(3)
					prev_rod_pos = rod_pos
					rod_pos = bd%interface_volume_element*((/i,j,k/))+rod_start
			
					!evaluate p
					p=0
			
					do m=1, bd%NO
						diff = fr%molecule(m)%o%position - rod_pos
						!pbc correction              
						do s=1,2
							if (diff(s) > bd%box_dimensions(s)/2.0) then
								diff(s) = diff(s) - bd%box_dimensions(s)
							else if (diff(s) < -bd%box_dimensions(s)/2.0) then
								diff(s) = bd%box_dimensions(s) + diff(s)
							end if
						end do
	
						r = norm2(diff)
						if (r<=3*E) then !after 3 sigma cut off ... the value would be too small to include and it is expesnive to calculate exp...
							p = p + exp(-r**2/(2*E**2))/((2*pi*E**2)**1.5)
						end if  
					end do
			
					!shaping the p
					prev_pdiff = pdiff
					pdiff = abs(bd%interface_density*waterDensity - p)
			
					!if the positive derivative is found ams the function has value lower than tollerance which is a number close to zero -> we have found the interface
					if(prev_pdiff < tollerance) then
						if(pdiff > prev_pdiff) then
							found_bottom_interface = .true.
							points(i,j,2) = prev_rod_pos(3)
							prev_index(i,j,2) = k - 1
							exit
						end if
					end if
				end do
				
				rod_pos = rod_start
				pdiff = bd%interface_density*waterDensity
	
				!poking rod from top
				do k=min(n_points(3), prev_index(i,j,1)+bd%INTERFACE_PUSHBACK),0,-1
					prev_rod_pos = rod_pos
					rod_pos = bd%interface_volume_element*((/i,j,k/))+rod_start
			
					!evaluate p
					p=0
			
					do m=1, bd%NO
						diff = fr%molecule(m)%o%position - rod_pos
						!pbc correction              
						do s=1,2
							if (diff(s) > bd%box_dimensions(s)/2.0) then
								diff(s) = diff(s) - bd%box_dimensions(s)
							else if (diff(s) < -bd%box_dimensions(s)/2.0) then
								diff(s) = bd%box_dimensions(s) + diff(s)
							end if
						end do
	
						r = norm2(diff)
						if (r<=3*E) then !after 3 sigma cut off ... the value would be too small to include and it is expesnive to calculate exp...
							p = p + exp(-r**2/(2*E**2))/((2*pi*E**2)**1.5)
						end if  
					end do
			
					prev_pdiff = pdiff
					pdiff = abs(bd%interface_density*waterDensity - p)
	
					!if the positive derivative is found ams the function has value lower than tollerance which is a number close to zero -> we have found the interface
					if(prev_pdiff < tollerance) then
						if(pdiff > prev_pdiff) then
							found_uper_interface = .true.
							points(i,j,1) = prev_rod_pos(3)
							prev_index(i,j,1) = k + 1
							exit
						end if
					end if
				end do
				
				!checks if program found both interfaces
		
				if((found_uper_interface .and. found_bottom_interface) == .false.) then
					if((found_uper_interface .or. found_bottom_interface) == .false.) then
						errmsg = " program did not find neither the lower nor the upper interface"
					else if(found_bottom_interface == .false.) then
						errmsg = " program did not find the lower interface"
					else if(found_uper_interface == .false.) then
						errmsg = " program did not find the upper interface"
					end if
			
					print "(A,i10,A)", "in step n: ", step+1, errmsg
			
					stop
				end if 
				
			end do
		
		end do
		!$OMP END DO
		!$OMP END PARALLEL
	end if
	
	call frame_density_function()
	!instead of skipping frames calculate density file with all of them
	do k = 1, min(bd%INTERFACE_SKIP-1, bd%NSTEP-step)
		print "(a,a,i,$)", char(13), "STEP: ", step+k
		call fr%read_frame()
		call frame_density_function()
	end do
	
	!printing out each step in interface.xyz
	if(vmdout) then 
		call vmd_out()
	end if
	
	if(bd%cancel_interface_calculation == .false.) then
		call bin_out()
	end if
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end do

if(vmdout) then 
	close(f_interface)
end if

close(f_interfacebin)

!normalize the density function and write it into files
call evaluate_density_function()

deallocate(points)
deallocate(density_function)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS AND SUBROUTINES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains

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
					print*, "ERROR: Interface program is not supported for hydroxyl trajectories"
					stop
				end if
				i = i + 1
			case('g')
				nocheck = .true.
				i = i + 1
			case('v')
				vmdout = .true.
				i = i + 1
			case default
				print*, "skipping invalid option ", trim(op)
				i = i + 1
		end select

		! all the switches were evaluated
		if(i > iargc()) exit
	end do

	print*, ""

	! check if all the mandatory options were selected
	if(.not. associated(fr)) then
		print*, "input file(s) must be specified, see help -h"
		stop
	end if

end subroutine evaluate_switches

subroutine printHelp()

	print "(a)",		""
	print "(a)",		"__________________________Help dialog of Interface_____________________________"
	print "(a,tr8,a)",	""                          
	print "(a,tr8,a)",	"_____________________________Mandatory options_________________________________"
	print "(a,tr8,a)",	""
	print "(a)",		"-I(input) <file1> *<file2>"
	print "(tr4,a)",		"specifies the input files. Options are:"
	print "(tr8,a)",			"*.xyz (positions), *.xyz (velocities)"
	print "(tr8,a)",			"*.gro (positions+velocities)"
	print "(tr8,a)",			"*.trr, *.gro (one frame - can be only positions)"
	print "(a)",		""
	print "(a)",		"___________________________Optional useful options_____________________________"
	print "(a)",		""
	print "(a)",	    "-G(good)"
	print "(tr4,a)",	    "skips checking of the input files"
	print "(a)",	    ""
	print "(a)",		"-H(help)"
	print "(tr4,a)",		"prints this help dialog"
	print "(a)",		""
	print "(a)",		"__________________________Optional debugging options___________________________"
	print "(a)",		""
	print "(a)",		"WARNING:"
	print "(tr4,a)",		"These options will create HUGE ASCII files."
	print "(tr4,a)",		"These options are kept in the code mainly for historical reasons/debugging."
	print "(a)",		""
	print "(a)",		"-V(vmdout)"
	print "(tr4,a)",		"ASCII dump of interface points in .xyz format"
	print "(a)",		""
	print "(a)", "for more details see documentation"

end subroutine printHelp

function abs_pbc_check(vector,box_dimensions)
	real*8, dimension(3) :: abs_pbc_check
	real*8, dimension(3), intent(IN) :: vector
	real*8, dimension(3), intent(IN) :: box_dimensions
	integer :: i
	
	abs_pbc_check = vector
	
	do i=1,3
		if (abs_pbc_check(i) < 0) then
			abs_pbc_check(i) = -abs_pbc_check(i)
		end if
		if (abs_pbc_check(i) > box_dimensions(i)/2) then
			abs_pbc_check(i) = box_dimensions(i) - abs_pbc_check(i)
		end if
	end do
	
end function abs_pbc_check


subroutine poke_rod() !private variables for parralelization should be: k,prev_rod_pos, rod_pos, p, pdiff, prev_pdiff , m, diff 
	!poking rod from bottom
	!k = maxval((/int(0, 8),prev_index(i,j,2)-20/)) <- seems weird, but max() does not work??? error #6414 this PARAMETER constant name is invalid in this context
	do k=maxval((/int(0, 8),prev_index(i,j,2)-20/)) ,n_points(3)
		prev_rod_pos = rod_pos
		rod_pos = bd%interface_volume_element*((/i,j,k/))+rod_start
		
		!evaluate p
		p=0
		
		do m=1, bd%NO
			diff = fr%molecule(m)%o%position - rod_pos
			!pbc correction              
			diff = abs_pbc_check(diff, bd%box_dimensions)
					
			r = norm2(diff)
			if (r<=3*E) then !after 3 sigma cut off ... the value would be too small to include and it is expesnive to calculate exp...
				p = p + exp(-r**2/(2*E**2))/((2*pi*E**2)**1.5)
			end if  
		end do
		
		!shaping the p
		prev_pdiff = pdiff
		pdiff = abs(bd%interface_density*waterDensity - p)
		
		!if the positive derivative is found ams the function has value lower than tollerance which is a number close to zero -> we have found the interface
		if(pdiff < tollerance) then
			if(pdiff > prev_pdiff) then
				found_bottom_interface = .true.
				points(i,j,2) = prev_rod_pos(3)
				prev_index(i,j,2) = k
				exit
			end if
		end if
	end do
			
	!poking rod from top
	!k=min((/n_points(3), prev_index(i,j,1)+20/)) <- seems weird, but max() does not work??? error #6414 this PARAMETER constant name is invalid in this context
	do k=minval((/n_points(3), prev_index(i,j,1)+20/)),0,-1
		prev_rod_pos = rod_pos
		rod_pos = bd%interface_volume_element*((/i,j,k/))+rod_start
		
		!evaluate p
		p=0
		
		do m=1, bd%NO
			diff = fr%molecule(m)%o%position - rod_pos
			!pbc correction              
			diff = abs_pbc_check(diff, bd%box_dimensions)
					
			r = norm2(diff)
			if (r<=3*E) then !after 3 sigma cut off ... the value would be too small to include and it is expesnive to calculate exp...
				p = p + exp(-r**2/(2*E**2))/((2*pi*E**2)**1.5)
			end if  
		end do
		
		prev_pdiff = pdiff
		pdiff = abs(bd%interface_density*waterDensity - p)
		
		!if the positive derivative is found ams the function has value lower than tollerance which is a number close to zero -> we have found the interface
		if(pdiff < tollerance) then
			if(pdiff > prev_pdiff) then
				found_uper_interface = .true.
				points(i,j,1) = prev_rod_pos(3)
				prev_index(i,j,1) = k
				exit
			end if
		end if
	end do
end subroutine poke_rod

subroutine d_check()
	!checks if program found both interfaces
	
	if((found_uper_interface .and. found_bottom_interface) == .false.) then
		if((found_uper_interface .or. found_bottom_interface) == .false.) then
			errmsg = " program did not find neither the lower nor the upper interface"
		else if(found_bottom_interface == .false.) then
			errmsg = " program did not find the lower interface"
		else if(found_uper_interface == .false.) then
			errmsg = " program did not find the upper interface"
		end if
		
		print "(A,i10,A)", "in step n: ", step, errmsg
		
		stop
	end if 
	
end subroutine d_check

subroutine vmd_out()
	rod_start = -bd%box_dimensions/2.0 + bd%interface_volume_element/2.0
	do i = 0, bd%INTERFACE_SKIP
		s=step+i
		write(f_interface,*) size(points)
		write(f_interface,*) "step =", s
		do j = 0, n_points(1)
			do k = 0, n_points(2)
				do m = 1, 2
					write(f_interface,*) "p", (/bd%interface_volume_element(1)*(j)+rod_start(1), bd%interface_volume_element(2)*(k)+rod_start(2), points(j,k,m)/)
				end do
			end do
		end do
	end do
end subroutine vmd_out

subroutine bin_out()
	!write header at first!
	if(bin_header_done == .FALSE.) then
		!number of points in each step...
		write(f_interfacebin) int(size(points),8)
			!x,y grid
			do j = 0, n_points(1)
				do k = 0, n_points(2)
					write(f_interfacebin) bd%interface_volume_element(1)*(j)+rod_start(1), bd%interface_volume_element(2)*(k)+rod_start(2)
				end do
			end do
		bin_header_done = .TRUE.
	end if
	
	!continue with Z values bottom(1), top(2)
	
	do j = 0, n_points(1)
		do k = 0, n_points(2)
			do m = 1, 2
				write(f_interfacebin) points(j,k,m)
			end do
		end do
	end do

end subroutine bin_out

!density file subs and fs

function pbc_check(vector,box_dimensions)

	real*8, dimension(3) :: pbc_check
	real*8, dimension(3), intent(IN) :: vector
	real*8, dimension(3), intent(IN) :: box_dimensions
	integer :: i
	
	pbc_check = vector
	
	do i=1,3
		if (pbc_check(i) > box_dimensions(i)/2) then
			pbc_check(i) = pbc_check(i) - box_dimensions(i)
		else if (pbc_check(i) < -box_dimensions(i)/2) then
			pbc_check(i) = box_dimensions(i) + pbc_check(i)
		end if
	end do

end function pbc_check

subroutine frame_density_function
	integer :: x,y,x1,y1 !the position in the grid
	real*8 :: uz, dz !interface z value at position
	real*8 :: xgrid, ygrid, xt, yu
	real*8, dimension(3) :: diff

	!$OMP PARALLEL DEFAULT(SHARED)
	!$OMP DO PRIVATE(m, diff, xgrid, ygrid, xt, yu, x, y, x1, y1, uz, dz) REDUCTION(+:density_function)
	do m=1,bd%NO
		diff = pbc_check(fr%molecule(m)%o%position, bd%box_dimensions)
		
		!find where the molecule belongs
		!X
		if( (diff(1) <= rod_start(1)) .or. (diff(1) > rod_start(1)+(n_points(1)+1)*bd%interface_volume_element(1)) ) then
			x = -1
		else
			x = ((diff(1)-rod_start(1)) / bd%interface_volume_element(1))
			x1 = x + 1
		end if
		!Y
		if( (diff(2) <= rod_start(2)) .or. (diff(2) > rod_start(2)+(n_points(2)+1)*bd%interface_volume_element(2)) ) then
			y = -1	
		else
			y = ((diff(2)-rod_start(2)) / bd%interface_volume_element(2))
			y1 = y + 1
		end if

		!now we know where on the grid the molecule lives...
		!time to get the interface height on the molecules position...
		
		xgrid = rod_start(1) + x*bd%interface_volume_element(1)
		ygrid = rod_start(2) + y*bd%interface_volume_element(2)

		!deal with the edge cases
		if(x < 0) then
			x=n_points(1)
			x1=0
		end if
		if(y < 0) then
			y=n_points(2)
			y1=0
		end if
		if(x1 > n_points(1)) then
			x1 = n_points(1)
			x = x1 - 1
		end if
		if(y1 > n_points(2)) then
			y1 = n_points(2)
			y = y1 - 1
		end if
		
		xt = (diff(1)-xgrid)/(bd%interface_volume_element(1))
		yu = (diff(2)-ygrid)/(bd%interface_volume_element(2))
		
		!interface z
		uz = (1-xt)*(1-yu)*points(x,y,1) + xt*(1-yu)*points(x1,y,1) + (1-xt)*yu*points(x,y1,1) + xt*yu*points(x1,y1,1) 
		dz = (1-xt)*(1-yu)*points(x,y,2) + xt*(1-yu)*points(x1,y,2) + (1-xt)*yu*points(x,y1,2) + xt*yu*points(x1,y1,2) 
		!distance from the interface
		uz = uz - diff(3)
		dz = diff(3) - dz
		
		! evaluating density function
		do i=1,nop
			if (uz >= (density_axis(i) - bd%density_tol) .and. (uz < density_axis(i)+bd%density_tol)) then
				density_function(i,1) = density_function(i,1) + 1
				density_function(i,2) = density_function(i,2) + 1
			end if
			if (dz >= (density_axis(i) - bd%density_tol) .and. (dz < density_axis(i) + bd%density_tol)) then
				density_function(i,1) = density_function(i,1) + 1
				density_function(i,3) = density_function(i,3) + 1
			end if
		end do
	end do
	!$OMP END DO
	!$OMP END PARALLEL

end subroutine frame_density_function

subroutine evaluate_density_function()
	integer ::                          f_density_vs_r
	
	open (newunit = f_density_vs_r, FILE="density.dat", RECL = 140, iostat = ierr)
	if(ierr .ne. 0) then
		print*, "ERROR: unable to create file called 'density.dat'" 
		stop
	end if
	
	write(f_density_vs_r,*) "#z-coordinate   av-density   up-density   down-density"
	
	do i=1,nop
		!print*, "dist: ", dist(i)
		!reporting number density in N/nm^3
		density_function(i,1) = ( density_function(i,1) / (bd%NSTEP*bd%density_tol*bd%box_dimensions(1)*bd%box_dimensions(2)*4) ) * 1000
		density_function(i,2) = ( density_function(i,2) / (bd%NSTEP*bd%density_tol*bd%box_dimensions(1)*bd%box_dimensions(2)*2) ) * 1000
		density_function(i,3) = ( density_function(i,3) / (bd%NSTEP*bd%density_tol*bd%box_dimensions(1)*bd%box_dimensions(2)*2) ) * 1000
		write(f_density_vs_r,*) density_axis(i) , density_function(i,1), density_function(i,2), density_function(i,3)
	end do
	
	close(f_density_vs_r)
	
end subroutine evaluate_density_function

subroutine read_interface_frame()
	real*8 :: placeholder
	integer*8 :: npoint
	
	if( (mod(step,bd%INTERFACE_SKIP) == 0) ) then
		if(bin_header_done == .false.) then
			read(f_interfacebin, iostat = ierr), npoint
			if(ierr .ne. 0) then
				print*, "ERROR: 'interface.bin' does not have enough data inside"
				stop
			end if
			
			if(npoint .ne. size(points)) then
				print*, "ERROR: 'interface.bin' does not match the current setup or is corrupted"
				stop
			end if

			!reading x,y values of points...
			do i = 1, size(points)/2
				read(f_interfacebin, iostat = ierr) placeholder, placeholder
				if(ierr .ne. 0) then
					print*, "ERROR: 'interface.bin' does not have enough data inside"
					stop
				end if
			end do
			bin_header_done = .true.
		end if
		
		!reading the z values of the points...
			
		do i = 0, n_points(1)
			do j = 0, n_points(2)
				read(f_interfacebin, iostat = ierr) points(i, j, 1), points(i, j, 2)
				if(ierr .ne. 0) then
					print*, "ERROR: 'interface.bin' does not have enough data inside"
					stop
				end if
			end do
		end do
	end if  

end subroutine read_interface_frame

end program interf
