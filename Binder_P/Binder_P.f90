program Binder_P
use bdata
use SFG_TYPES
use SFG_FRAMES
use switches

implicit none

type(boxdata) ::                                            bd

class(frame_reader), pointer ::                             fr
type(trr_frame_reader), allocatable, target ::              trr
type(gro_frame_reader), allocatable, target ::              gro
type(xyz_frame_reader), allocatable, target ::              xyz

real*8, dimension(:,:,:), allocatable ::                    points
!holds up whole frame of interface.xyz
!dimension(numberofpoints/2, 2 [upper,lower], 3 [x,y,z])

real*8, dimension(3) ::                                     diff
!holds calculations - xdiff,ydiff,zdiff

real*8, dimension(:), allocatable ::                        oih,&
															ojh,&
															ojoi
!oih, ojh oioj position vectors...

real*8 ::                                                   r,&
															lzdiff,&
															lznear,&
															rmin,&
															zcenter,&
															ojhdist,&
															oihdist,&
															oiojdist,&
															cosoho,&
															cosoonorm,&
															cosref,&
															distbinwidth,&
															anglebinwidth,&
															bottomrmin,&
															toprmin,&
															dipolecos
!r, lzdiff, lznear, rmin used to hold calculations
!oih ojh oo dists - for evaluation of hbonds
!cosoho - cos of angle between O-H..O 
!cosoonorm - cos of O->O and surface normal (pointing from the interface to the bulk)
!cosref - just storing value of cos(140 deg)
!dist start-end, angle start-end - boundaries for 3d histogram (O-O dist & cos O->O vs surface normal) distribution...
!distbinwidth anglebinwidth - bin widths for the 3D histogram
!bottom & top rmin - shortest distance of a molecule to a bottom & top interface point
!bd%dipole_rstart, bd%dipole_rend, bd%dipole_rbinwidth ... parameters for the dipole moment plot...
!dipolecos - value of the cosine dipole moment versus interface normal...

integer*8, dimension(0:15) :: layercount,steplayercount

real*8, dimension(0:15,0:15) :: cumulativedonorcount, cumulativeacceptorcount

real*8, dimension(:,:), allocatable :: dipoleanglevsr
!size (2, ...) holds dipole angle versus r

integer*8, dimension(:,:,:,:), allocatable ::                   histogram
!array to hold counts for the 3d histogram (type, inter/intra, binder ,0:bd%hbhist_distdiv-1, 0:bd%hbhist_anglediv-1)
!
!inter/intra:
!1 intra
!2 inter
!
!binder:
!0-15 same meaning as binder...
!
!0:distsiv...:
!size in X axis
!
!0:bd%hbhist_anglediv...:
!size in Y axis


integer*8, dimension(:,:), allocatable ::                   donorCountMatrix,&
															dipoleanglevsrcounts

! donorCountMatrix - notation donor-acceptor - for example dL0-dL1 means dL0 donor dL1 acceptor
!matrix looks like this:
! dL0-dL0     dL0-dL1     dL0-dL2     ...     dL0-uL3
! dL1-dL0     dL1-dL1       ...       ...     dL1-uL3
!   .                        .                   .
!   .                                  .         .
!   .                                            .
! uL3dL0       ...          ...       ...     uL3-uL3

!dipoleanglevsrcounts - used to average the plot based on the number of contributions to each bin

integer*1, dimension(:), allocatable :: binder

integer*1, dimension(:), allocatable ::                     HBONDFileUnits,&
															histogramFileUnits
!binder stores layer affiliation of each water molecule
! 0 uL0     4 dL0
! 1 uL1     5 dL1
! 2 uL2     6 dL2
! 3 uL3     7 dL3

integer*8 ::                                                step,&
															i,&
															j,&
															k,&
															m,&
															n,&
															Npoint,&
															r_index
!step - iterator through steps
!i,j,m,k -general iterators
!npoint - number of points in interface.xyz file
!bd%hbhist_distdiv bd%hbhist_anglediv - number of divisions in range dist or angle (end - start)

!uL0 ... count of oxygens belonging to each layer - mk array

integer ::                                                  ierr,&
															f_interfacebin,&
															f_newbinder,&
															f_hbonds,&
															f_dipole

integer*1 ::                                                upordown,&
															outbyte
!upordown - oxygen belongs to bottom layer --- 0
!upordown - oxygen belongs to upper layer --- 1

!outbyte - byte to be printed into file
!0 1 2 3 - upper layer 0 1 2 3
!4 5 6 7 - bottom layer 0 1 2 3
!easy to compare with (upperorlower*4 + layernumber) in later steps

character(256)                                              txt

character(len=60), allocatable ::                           arg,&
															arg1

logical ::                                                  hbhist = .false.,&
															dipole = .false.,&
															bin_header_done = .FALSE.,&
															nocheck = .false.,&
															layerFound = .false.

!technical debt variables
real*8, dimension(3) ::										rod_start

integer*8, dimension(3) ::									n_points

integer*8 ::												points_index00,&
															points_index01,&
															points_index10,&
															points_index11,&
															x,&
															x1,&
															y,&
															y1

real*8 ::													xt,&
															yu,&
															xgrid,&
															ygrid

!##########################################EVALUATE PROGRAM SWITCHES###########################################

call evaluate_switches()
	
!##############################################################################################################

!! READING BOXDATA

call bd%read_boxdata()

!!!END OF READING BOXDATA

!!check layers_limits...

do i=2,bd%get_layersCount()-1
	if(bd%layers_limits(i-1) >= bd%layers_limits(i)) then
		print*, "ERROR: $LAYERS_LIMITS are not in ascending order"
		print*, "$LAYERS_LIMITS :"
		print*, bd%LAYERS_LIMITS
		stop
	end if
end do

if( (hbhist == .true.) .or. (dipole == .true.) ) then
	print*, "Recap of parameters:"
	print*, ""
end if

if(hbhist == .true.) then
	print*, "-hbhist recap:"
	print*, "range of O-O distance:"
	print*, "from ", bd%hbhist_diststart , "to ", bd%hbhist_distend
	print*, "number of bins in this range: ", bd%hbhist_distdiv
	print*, "range of surface normal-O-O cos(angle):"
	print*, "from ",bd%hbhist_anglestart , "to ", bd%hbhist_angleend
	print*, "number of bins in this range: ", bd%hbhist_anglediv
	print*, ""
end if

if(dipole == .true.) then
	print*, "-dipole recap:"
	print*, "range of distance from the instantaneous surface:"
	print*, "from ", bd%dipole_rstart , "to ", bd%dipole_rend
	print*, "width of bins: ", bd%dipole_rbinwidth 
	print*, ""
end if

!!opening xyz files
!allocate(xyz)
!fr => xyz
!call fr%open_file(30, "pos_rebuilt.xyz", "vel_mol_rebuilt.xyz")
allocate(binder(bd%NO))


if(.not. nocheck) then
	print*, "Checking the contents of trajectory file(s)"
	print*, ""
	do i = 1, bd%NSTEP
		call fr%skip_frame()
	end do
	call fr%rewind_file()
	print*, "contents of trajectory files - OK"
	print*, ""
end if

!checking if interface.bin has all the points needed...
print*, "Checking 'interface.bin' file"

open(newunit = f_interfacebin, file = "interface.bin", form = "unformatted", access = 'stream', convert = 'big_endian', status = 'old', iostat = ierr)
if(ierr .ne. 0) then
	print*, "ERROR: unable to open file 'interface.bin'" 
	stop
end if

read(f_interfacebin, iostat = ierr) npoint
if(ierr .ne. 0) then
	print*, "ERROR: 'interface.bin' does not have enough data inside"
	stop
end if

allocate(points(Npoint/2,2,3))

rewind(f_interfacebin)

!checking if interfacebin has all the frames needed

do step = 0, bd%NSTEP-1
	if( (mod(step,bd%INTERFACE_SKIP) == 0) ) then
		if(bin_header_done == .false.) then
			read(f_interfacebin, iostat = ierr), npoint
			if(ierr .ne. 0) then
				print*, "ERROR: 'interface.bin' does not have enough data inside"
				stop
			end if
			
			!reading x,y values of points...
			do i = 1, npoint/2
				read(f_interfacebin, iostat = ierr) points(i, 1, 1), points(i, 1, 2)
				if(ierr .ne. 0) then
					print*, "ERROR: 'interface.bin' does not have enough data inside"
					stop
				end if
				points(i, 2, :) = points(i, 1, :)
			end do
			bin_header_done = .true.
		end if
		
		!reading the z values of the points...
			
		do i = 1, npoint/2
			read(f_interfacebin, iostat = ierr) points(i, 1, 3), points(i, 2, 3)
			if(ierr .ne. 0) then
				print*, "ERROR: 'interface.bin' does not have enough data inside"
				stop
			end if
		end do
	end if    
end do

!resetting
bin_header_done = .false.
points = 0
print*, "file 'interface.bin' - OK"
print*, ""

close(f_interfacebin)

open(newunit = f_interfacebin, file = "interface.bin", form = "unformatted", access = 'stream', convert = 'big_endian', status = 'old', iostat = ierr)
if(ierr .ne. 0) then
	print*, "ERROR: unable to open file 'interface.bin'" 
	stop
end if

!NEW BINDER
open(newunit = f_newbinder, FILE="binder.bin", form = "unformatted", access="stream", convert = 'big_endian', iostat = ierr)
if(ierr .ne. 0) then
	print*, "ERROR: unable to create file 'binder.bin'" 
	stop
end if

do step = 0, bd%NSTEP-1
	print "(a,a,i,$)", char(13), "STEP: ", step+1
	!reading frame
	call fr%read_frame()
	
	!estimating center of the water slab
	if (step == 0) then
		print*, ""
		print*, "estimating center of the water slab..."
		zcenter = 0
		do i = 1, bd%NO
			zcenter = zcenter + fr%molecule(i)%o%position(3)
		end do
		zcenter = zcenter/bd%NO !searching center of water slab as center of mass of all the oxygens
		print*, "estimated center of the water slab is at Z = ", zcenter
	end if
	
	! reading interface points
	if( (mod(step,bd%INTERFACE_SKIP) == 0) ) then
		if(bin_header_done == .false.) then
			!first interface read ever
			read(f_interfacebin), npoint
			
			n_points = NINT(bd%box_dimensions/bd%interface_volume_element)

			bd%interface_volume_element = bd%box_dimensions/n_points

			rod_start = -bd%box_dimensions(:)/2.0 + bd%interface_volume_element(:)/2.0
			
			!reading x,y values of points...
			do i = 1, npoint/2
				read(f_interfacebin) points(i, 1, 1), points(i, 1, 2)
				points(i, 2, :) = points(i, 1, :)
			end do
			bin_header_done = .true.
		end if
		
		!reading the z values of the points...
			
		do i = 1, npoint/2
			read(f_interfacebin) points(i, 1, 3), points(i, 2, 3)
		end do
		
	end if
	
	!prepare for dipole plots...
	if( (step == 0) .and. dipole) then
			allocate(dipoleanglevsr(2, 0:int((bd%dipole_rend-bd%dipole_rstart)/bd%dipole_rbinwidth)-1))
			dipoleanglevsr = 0
			allocate(dipoleanglevsrcounts(2, 0:int((bd%dipole_rend-bd%dipole_rstart)/bd%dipole_rbinwidth)-1))
			dipoleanglevsrcounts = 0
	end if
	
	! check distance of each  molecules from the interface
	!$OMP PARALLEL DEFAULT(SHARED)
	!$OMP DO PRIVATE(r_index, rmin, bottomrmin, toprmin, diff, r, i, j, m, lznear, dipolecos, lzdiff, upordown, layerFound, xt, yu, points_index00, points_index01, points_index10, points_index11, xgrid, ygrid, x, x1, y, y1 ), REDUCTION(+:dipoleanglevsr, dipoleanglevsrcounts)
	do m=1, bd%NO
		
		! new density approach
		! the whole purpose of this block is to feed the bottomrmin and toprmin variables
		! todo here will be a technical debt, and I dont care since this whole program is a mess...
		! the bottomrzmin and toprzmin makes little to no sense... so it will be set to -inf...
		diff = pbc_check(fr%molecule(m)%o%position, bd%box_dimensions)
	
		!find where the molecule belongs
		!X
		if( (diff(1) <= rod_start(1)) .or. (diff(1) > rod_start(1)+(n_points(1))*bd%interface_volume_element(1)) ) then
			x = -1
		else
			x = ((diff(1)-rod_start(1)) / bd%interface_volume_element(1))
			x1 = x + 1
		end if
		!Y
		if( (diff(2) <= rod_start(2)) .or. (diff(2) > rod_start(2)+(n_points(2))*bd%interface_volume_element(2)) ) then
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
			x=n_points(1)-1
			x1=0
		end if
		if(y < 0) then
			y=n_points(2)-1
			y1=0
		end if
		if(x1 > n_points(1)-1) then
			x1 = n_points(1)-1
			x = x1 - 1
		end if
		if(y1 > n_points(2)-1) then
			y1 = n_points(2)-1
			y = y1 - 1
		end if

		points_index00 = x*n_points(2) + y + 1
		points_index01 = x*n_points(2) + y1 + 1
		points_index10 = x1*n_points(2) + y + 1
		points_index11 = x1*n_points(2) + y1 + 1
		
		xt = (diff(1)-xgrid)/(bd%interface_volume_element(1))
		yu = (diff(2)-ygrid)/(bd%interface_volume_element(2))
		
		!interface z
		toprmin = (1-xt)*(1-yu)*points(points_index00,1,3) + xt*(1-yu)*points(points_index10,1,3) + (1-xt)*yu*points(points_index01,1,3) + xt*yu*points(points_index11,1,3) 
		bottomrmin = (1-xt)*(1-yu)*points(points_index00,2,3) + xt*(1-yu)*points(points_index10,2,3) + (1-xt)*yu*points(points_index01,2,3) + xt*yu*points(points_index11,2,3) 
		!distance from the interface
		toprmin = toprmin - diff(3)
		bottomrmin = diff(3) - bottomrmin

		!each molecule will add a contribution to dipole angle vs rmin...
		
		!cos of dipole vs [0,0,1] calculation...
		if(dipole) then
			!!!
			diff = pbc_check(fr%molecule(m)%o%position(:)-fr%molecule(m)%h1%position(:), bd%box_dimensions)
			diff = diff + pbc_check(fr%molecule(m)%o%position(:)-fr%molecule(m)%h2%position(:), bd%box_dimensions)
			diff = diff/norm2(diff)
			!diff is now unit vector in the direction of water dipole...
			dipolecos = dot_product(diff, (/0,0,1/))
			!the dipole was flipped...
			dipolecos = -dipolecos
			!!!
		
			!from bottom interface
		
			!vs r from interface...
			r_index = int((bottomrmin-bd%dipole_rstart)/bd%dipole_rbinwidth)
			!record only if in the boundaries of the array...
			if( (r_index <= int((bd%dipole_rend-bd%dipole_rstart)/bd%dipole_rbinwidth)-1) .and. (r_index >= 0) ) then 
				dipoleanglevsr(1,r_index) = dipoleanglevsr(1,r_index) + dipolecos
		
				dipoleanglevsrcounts(1,r_index) = dipoleanglevsrcounts(1,r_index) + 1
			end if

			!from top interface
		
			!vs r from interface...
			r_index = int((toprmin-bd%dipole_rstart)/bd%dipole_rbinwidth)
			!record only if in the boundaries of the array...
			if( (r_index <= int((bd%dipole_rend-bd%dipole_rstart)/bd%dipole_rbinwidth)-1) .and. (r_index >= 0) ) then 
				dipoleanglevsr(2,r_index) = dipoleanglevsr(2,r_index) + dipolecos
		
				dipoleanglevsrcounts(2,r_index) = dipoleanglevsrcounts(2,r_index) + 1
			end if

		end if
		
		!binder label assignment
		if(fr%molecule(m)%o%position(3) >= zcenter) then
			upordown = 1
			rmin = toprmin
		else 
			upordown = 0
			rmin = bottomrmin
		end if
		
		layerFound = .false.
		! is layer 0?
		if ( rmin <= bd%layers_limits(1) ) then
			!print*, "molecule ", m, "belongs to L 0 ", upordown
			binder(m) = 0+upordown*8
			layerFound = .true.
		end if
		
		! is it layer < layers_count?
		if( .not. layerFound ) then
			do i=2, bd%get_layersCount()-1
				!layer i
				if (rmin > bd%layers_limits(i-1) .and. rmin <= bd%layers_limits(i) ) then
					binder(m) = (i-1)+upordown*8
					layerFound = .true.
					exit
				end if
			end do
		end if
		
		! it must be layer == layersCount
		if( .not. layerFound) then
			binder(m) = (bd%get_layersCount()-1)+upordown*8
		end if
		
	end do
	!$OMP END DO
	!$OMP END PARALLEL
		
	!write a frame into binder.bin
	do m=1, bd%NO
		
		if(mod(m,2)==1) then
			outbyte = 0
			!I need to make 4 MSB of outbyte to be equal to binders 4 LSB (basically value 0-15)
			outbyte = IOR(outbyte, ISHFT(binder(m),4))
		else 
			! need to make 4 LSB of outbyte to be equal to binders 4LSB
			outbyte = IOR(outbyte, IAND(binder(m), Z'0F'))
		end if
		
		!write byte every 2 molecules
		if(mod(m,2)==0) then
			write(f_newbinder) outbyte
			!print*, RSHIFT(outbyte, 4), (outbyte .AND. X'0F') !THIS IS HOW TO READ THE FILE BYTE BY BYTE
		end if
		
		!if odd number of molecules write the last byte with 0 in top nibble
		if(m == bd%NO) then
			if(mod(bd%NO,2)) then
				write(f_newbinder) outbyte
			end if    
		end if
	end do
	
	!binder is done for this frame
	
	!time to process the HBONDS
	!using i,j for water molecules
	
	!allocating memory for hbhists
	if(hbhist) then
		if( (mod(step,bd%HBHIST_SKIP) == 0) ) then
			if(step == 0) then
				distbinwidth = (bd%hbhist_distend - bd%hbhist_diststart)/bd%hbhist_distdiv
				anglebinwidth = (bd%hbhist_angleend - bd%hbhist_anglestart)/bd%hbhist_anglediv
			
				allocate(oih(3))
				allocate(ojh(3))
				allocate(ojoi(3))
				allocate(histogram(2,0:15,0:bd%hbhist_distdiv-1, 0:bd%hbhist_anglediv-1))
				histogram = 0
				allocate(donorCountMatrix(0:15,0:15))
				donorCountMatrix = 0
				cosref = 7d0*dacos(-1d0)/9d0 !7pi/9 = 140 deg
				cosref = dcos(cosref)
				layercount = 0
				call deployHistogramFiles()
				
				cumulativedonorcount = 0
				cumulativeacceptorcount = 0
			end if
		
			steplayercount = 0
			donorcountmatrix = 0
		
			!paralelize this loop...
	!$OMP PARALLEL DEFAULT(shared)
		!$OMP DO PRIVATE(j, ojoi, oiojdist, m, oih, ojh, oihdist, ojhdist, cosoho, cosoonorm) REDUCTION(+:donorcountmatrix, layercount, steplayercount)
		!firstprivate cosref
		!private j, ojoi, oiojdist, m, oih, ojh, oihdist, ojhdist, cosoho, cosoonorm
		!reduction donorcountmatrix, layercount, steplayercount
			do i=1, BD%NO
				layercount(binder(i)) = layercount(binder(i)) + 1
				steplayercount(binder(i)) = steplayercount(binder(i)) + 1
				do j=1, BD%NO
					if(i .NE. j) then
						!first check the oo distance ojh is oi<-oj in this case
						ojoi = fr%molecule(j)%o%position - fr%molecule(i)%o%position
						ojoi = pbc_check(ojoi, bd%box_dimensions)
						oiojdist = norm2(ojoi)
					
						if(oiojdist < 3.2) then !if oo dist < 3.2 then check for angle
							do m=1, 2 !for both O-H in water molecule
								if(mod(m,2) .EQ. 1) then
									oih = fr%molecule(i)%h1%position - fr%molecule(i)%o%position !oi->h1
									ojh = fr%molecule(i)%h1%position - fr%molecule(j)%o%position !oj->h1
									oih = pbc_check(oih, bd%box_dimensions)
									ojh = pbc_check(ojh, bd%box_dimensions)
								else
									oih = fr%molecule(i)%h2%position - fr%molecule(i)%o%position !oi->h2
									ojh = fr%molecule(i)%h2%position - fr%molecule(j)%o%position !oj->h2
									oih = pbc_check(oih, bd%box_dimensions)
									ojh = pbc_check(ojh, bd%box_dimensions)
								end if
							
								oihdist = norm2(oih)
								ojhdist = norm2(ojh)

								!now checking the angle...
								cosoho = dot_product(oih,ojh)/(oihdist*ojhdist)
								if(cosoho < cosref) then
									cosoonorm = dot_product(ojoi, surfaceNormal(binder(i)))/(oiojdist)
								
									!counting into donorAceeptorCountMatrix
									donorCountMatrix(binder(i),binder(j)) = donorCountMatrix(binder(i),binder(j)) + 1
									!counting for the histogram...
									if(oiojdist >= bd%hbhist_diststart) then !then the value belongs to the histogram
										!this is pretty ugly... fortran give me increment operator pls
									
										!hb donors...
										if(binder(i) == binder(j)) then !add 1 to intra layer histogram
											histogram(1,&
												binder(i),&
												int((oiojdist-bd%hbhist_diststart)/distbinwidth),&
												int((cosoonorm-bd%hbhist_anglestart)/anglebinwidth)&
											) = &
												histogram(1,&
													binder(i),&
													int((oiojdist-bd%hbhist_diststart)/distbinwidth),&
													int((cosoonorm-bd%hbhist_anglestart)/anglebinwidth)&
												) + 1
										else !else add 1 to inter binder(i) layer...
											histogram(2,&
												binder(i),&
												int((oiojdist-bd%hbhist_diststart)/distbinwidth),&
												int((cosoonorm-bd%hbhist_anglestart)/anglebinwidth)&
											) = &
												histogram(2,&
													binder(i),&
													int((oiojdist-bd%hbhist_diststart)/distbinwidth),&
													int((cosoonorm-bd%hbhist_anglestart)/anglebinwidth)&
												) + 1
										end if
									
									end if
								
								end if
					
							end do    
						end if
						
					end if
				end do
			end do
		!$OMP END DO
	!$OMP END PARALLEL
			do i = 0, 15
				do j = 0, 15
					if (steplayercount(i) .ne. 0) then
						cumulativedonorcount(i,j) = cumulativedonorcount(i,j) + real(donorcountmatrix(i,j),8)/steplayercount(i) 
						cumulativeacceptorcount(i,j) = cumulativeacceptorcount(i,j) + real(donorcountmatrix(j,i),8)/steplayercount(i)
					end if     
				end do
			end do
			
		end if
	end if	
end do

close(f_newbinder)

print*, ""

!these are simones total hbond numbers... average number of hbonds per molecule in a layer which can be even split into layerwise numbers...

if(hbhist) then
	open(newunit = f_hbonds, file = "hbonds.dat", iostat = ierr)
	if(ierr .ne. 0) then
		print*, "ERROR: unable to create file 'hbonds.dat'" 
		stop
	end if
	
	do i = 0, 1
		do j = 0, bd%get_layersCount()-1
			write(f_hbonds,"(A,TR3,I)") "LAYER:", j+8*i
			write(f_hbonds,"(A,F)") "average number of water molecules: ", real(layercount(j+8*i),8)/INT(bd%NSTEP/BD%HBHIST_SKIP,4)
			write(f_hbonds,"(A)") "layer totals:"
			write(f_hbonds,"(A,TR20,A,TR18,A)") " total", "donors", "acceptors"
			write(f_hbonds,*) sum(cumulativedonorcount(j+8*i,:))/bd%NSTEP + sum(cumulativeacceptorcount(j+8*i,:))/INT(bd%NSTEP/BD%HBHIST_SKIP,4),&
					sum(cumulativedonorcount(j+8*i,:))/INT(bd%NSTEP/BD%HBHIST_SKIP,4),&
					sum(cumulativeacceptorcount(j+8*i,:))/INT(bd%NSTEP/BD%HBHIST_SKIP,4)
			write(f_hbonds,"(A)") "intra layer totals:"
			write(f_hbonds,"(A,TR20,A,TR18,A)") " total", "donors", "acceptors"
			write(f_hbonds,*) cumulativedonorcount(j+8*i,j+8*i)/INT(bd%NSTEP/BD%HBHIST_SKIP,4) + cumulativeacceptorcount(j+8*i,j+8*i)/INT(bd%NSTEP/BD%HBHIST_SKIP,4),&
					cumulativedonorcount(j+8*i,j+8*i)/INT(bd%NSTEP/BD%HBHIST_SKIP,4),&
					cumulativeacceptorcount(j+8*i,j+8*i)/INT(bd%NSTEP/BD%HBHIST_SKIP,4)
			write(f_hbonds,"(A)") "inter layer totals:"
			write(f_hbonds,"(A,TR20,A,TR18,A)") " total", "donors", "acceptors"
			write(f_hbonds,*) (sum(cumulativedonorcount(j+8*i,:))-cumulativedonorcount(j+8*i,j+8*i))/INT(bd%NSTEP/BD%HBHIST_SKIP,4) + (sum(cumulativeacceptorcount(j+8*i,:))-cumulativeacceptorcount(j+8*i,j+8*i))/INT(bd%NSTEP/BD%HBHIST_SKIP,4),&
					(sum(cumulativedonorcount(j+8*i,:))-cumulativedonorcount(j+8*i,j+8*i))/INT(bd%NSTEP/BD%HBHIST_SKIP,4),&
					(sum(cumulativeacceptorcount(j+8*i,:))-cumulativeacceptorcount(j+8*i,j+8*i))/INT(bd%NSTEP/BD%HBHIST_SKIP,4)
			write(f_hbonds,"(A)") "respective layers:"
			write(f_hbonds,"(A,TR12,A,TR19,A,TR18,A)") "layer", "total", "donors", "acceptors"
		
			do m = 0, 1
				do n = 0, bd%get_layersCount()-1
					write(f_hbonds,"(I3,F,F,F)") n+m*8,&
								cumulativedonorcount(j+8*i,n+m*8)/INT(bd%NSTEP/BD%HBHIST_SKIP,4) + cumulativeacceptorcount(j+8*i,n+m*8)/INT(bd%NSTEP/BD%HBHIST_SKIP,4),&
								cumulativedonorcount(j+8*i,n+m*8)/INT(bd%NSTEP/BD%HBHIST_SKIP,4),&
								cumulativeacceptorcount(j+8*i,n+m*8)/INT(bd%NSTEP/BD%HBHIST_SKIP,4)
				end do
			end do
		write(f_hbonds,"(A)") "________________________________________________________________________________"
		end do    
	end do
	close(f_hbonds)

	!print normalized histograms into files...
	do i=0,1
		do j=0,bd%get_layersCount()-1 !loop through all the layers
	
			write(HistogramFileUnits(j+i*8),*) "#O-O distances  angles  inter  intra	tot" 
			do m=0, bd%hbhist_distdiv-1 !loop through all distdivs
				do n=0, bd%hbhist_anglediv-1 !loop through all angledivs
					write(HistogramFileUnits(j+i*8),*) bd%hbhist_diststart+(2*m+1)*distbinwidth/2,&
												bd%hbhist_anglestart+(2*n+1)*anglebinwidth/2,&
												real(histogram(1,j+i*8,m,n),8)/sum(histogram(1,j+i*8,:,:)),&
												real(histogram(2,j+i*8,m,n),8)/sum(histogram(2,j+i*8,:,:)),&
												real(histogram(1,j+i*8,m,n)+histogram(2,j+i*8,m,n),8)/(sum(histogram(1,j+i*8,:,:))+sum(histogram(2,j+i*8,:,:)))
				end do
			end do
		end do
	end do
	
	call terminateHistogramFiles()
end if

if(dipole) then
	!normalizing the dipoleanglevsr...
	do i = 0, (bd%dipole_rend-bd%dipole_rstart)/bd%dipole_rbinwidth - 1
		do j = 1,2
			dipoleanglevsr(j,i) = dipoleanglevsr(j,i)/dipoleanglevsrcounts(j,i)
		end do
	end do

	!printing the dipoleanglevsr to files...
	open(newunit = f_dipole, file = "dipole.dat", iostat = ierr, recl = 120)
	if(ierr .ne. 0) then
		print*, "ERROR: unable to create file 'dipole.dat'" 
		stop
	end if
	write(f_dipole,*) "#r    bottom cos(angle;r)     top cos(angle;r)"
	do i = 0, (bd%dipole_rend-bd%dipole_rstart)/bd%dipole_rbinwidth - 1
		write(f_dipole,*) bd%dipole_rstart+(2*i+1)*bd%dipole_rbinwidth/2, dipoleanglevsr(1,i), dipoleanglevsr(2,i)
	end do
	close(f_dipole)
	
end if

!close interface
close(f_interfacebin)

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
					print*, "ERROR: Binder program is not supported for hydroxyl trajectories"
					stop
				end if
				i = i + 1
			case('g')
				nocheck = .true.
				i = i + 1
			case('b')
				hbhist = .true.
				i = i + 1
			case('d')
				dipole = .true.
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
	print "(a)",		"___________________________Help dialog of Binder_______________________________"
	print "(a)",		""                          
	print "(a)",		"_____________________________Mandatory options_________________________________"
	print "(a)",		""
	print "(a)",		"-I(input) <file1> *<file2>"
	print "(tr4,a)",		"specifies the input files. Options are:"
	print "(tr8,a)",			"*.xyz (positions), *.xyz (velocities)"
	print "(tr8,a)",			"*.gro (positions+velocities)"
	print "(tr8,a)",			"*.trr, *.gro (one frame - can be only positions)"
	print "(a)",		""
	print "(a)",		"___________________________Optional useful options_____________________________"
	print "(a)",		""
	print "(a)",		"-B(bonds)"
	print "(tr4,a)",		"generate hydrogen bond distribution heatmaps for every layer"
	print "(a)",		""
	print "(a)",		"-D(dipole)"
	print "(tr4,a)",		"generate dipole orientation with respect to the instantaneous surface plots"
	print "(a)",		""
	print "(a)",		"-G(good)"
	print "(tr4,a)",		"skips checking of the input files"
	print "(a)",		""
	print "(a)",		"-H(help)"
	print "(tr4,a)",		"prints this help dialog"
	print "(a)",		""

end subroutine printHelp
	
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

subroutine deployHBONDFiles()
	allocate(HBONDFileUnits(0:15))
	HBONDFileUnits = 0
	
	!opening uL0 uL1 uL2 uL3 dL0 dL1 dL2 dL3 intra bond files
	!file unit is stored in HBONDFuleUnits on index 0 2 4 6 8 10 12 14
	do k = 0,1
		do i = 0, 3
			open(newunit=HBONDFileUnits(k*8+i*2), FILE=HBONDFileName(k*4+i,k*4+i), FORM='FORMATTED', iostat = ierr)
			if(ierr .ne. 0) then
				print"(a,a,a)", "ERROR: unable to create file '",trim(HBONDFileName(k*4+i,k*4+i)),"'"
				stop
			end if
			!write(HBONDFileUnits(k*8+i*2),*) "INTRA ", k*4+i, k*4+i
		end do
	end do
	!opening uL0 uL1 uL2 uL3 dL0 dL1 dL2 dL3 inter bond files
	!file unit is stored in HBONDFuleUnits on index 1 3 5 7 9 11 13 15
	do k = 0, 1
		do i = 0,3
			if(i == 3) then
				!getting layer indices 3,2 and 7,6
				open(newunit=HBONDFileUnits(k*8+i*2+1), FILE=HBONDFileName(k*4+i,k*4+i-1), FORM='FORMATTED', iostat = ierr)
				if(ierr .ne. 0) then
					print"(a,a,a)", "ERROR: unable to create file '",trim(HBONDFileName(k*4+i,k*4+i-1)),"'"
					stop
				end if
				!write(HBONDFileUnits(k*8+i*2+1),*) "INTER ", k*4+i, k*4+i-1
			else
				!getting indices 0,1 1,2 2,3 4,5 5,6 and 6,7
				!all the inter files are opened by now...
				open(newunit=HBONDFileUnits(k*8+i*2+1), FILE=HBONDFileName(k*4+i,k*4+i+1), FORM='FORMATTED', iostat = ierr)
				if(ierr .ne. 0) then
					print"(a,a,a)", "ERROR: unable to create file '",trim(HBONDFileName(k*4+i,k*4+i+1)),"'"
					stop
				end if
				!write(HBONDFileUnits(k*8+i*2+1),*) "INTER ", k*4+i, k*4+i+1
			end if
		end do
	end do
	
	print*, "donor files were deployed..."
end subroutine deployHBONDFiles

subroutine deployHistogramFiles()
	allocate(HistogramFileUnits(0:15))
	HistogramFileUnits = 0
	
	!opening histogram files for uL0 uL1 uL2 uL3 dL0 dL1 dL2 dL3 ... 
	do k = 0, 1
		do i = 0, bd%get_layersCount()-1
			open(newunit = HistogramFileUnits(k*8+i), FILE = HistogramFileName(k*8+i), FORM='FORMATTED', recl=360, iostat = ierr)
			if(ierr .ne. 0) then
				print"(a,a,a)", "ERROR: unable to create file '",trim(HistogramFileName(k*8+i)),"'"
				stop
			end if
		end do
	end do
	
	print*, "Histogram files were deployed..."
	
end subroutine deployHistogramFiles

subroutine terminateHBONDFiles()
	!closes all the HBOND files...
	do i = 0, 15
		close(HBONDFileUnits(i))
	end do

	deallocate(HBONDFileUnits)

	print*, "donor files were closed..."
end subroutine terminateHBONDFiles

subroutine terminateHistogramFiles()
	!closes all the histogram files...
	!todo if open close
	do i = 0, 15
		close(HistogramFileUnits(i))
	end do

	deallocate(HistogramFileUnits)
	print*, "Histogram files were closed..."
end subroutine terminateHistogramFiles

function getHBONDFileUnitIndex(binderi, binderj)
	!addition to deployHBONDFiles()
	!returns index of HBONDFileUnit based on binder value of two molecules
	integer, intent(IN) :: binderi, binderj
	integer :: getHBONDFileUnitIndex
	
	if(binderi == binderj) then
		getHBONDFileUnitIndex = binderi*2
	else
		getHBONDFileUnitIndex = (binderi*2)+1
	end if
	
end function getHBONDFileUnitIndex

function getHistogramFileUnitIndex(binderi)
	!addition to deployHBONDFiles()
	!returns index of HBONDFileUnit based on binder value of two molecules
	integer, intent(IN) :: binderi
	integer :: getHistogramFileUnitIndex
	
		getHistogramFileUnitIndex = binderi
	
end function getHistogramFileUnitIndex

function surfaceNormal(binderVal)
	real*8, dimension(3) :: surfaceNormal
	integer*1, intent(IN) :: binderVal
	
	if (binderVal < 8) then !bottom interface
		surfaceNormal = (/0,0,-1/)
	else !upper interface
		surfaceNormal = (/0,0,1/)
	end if 
	
end function surfaceNormal

function HBONDFileName(binderi, binderj)
	!addition to deployHBONDFiles()
	!returns name of the HBOND file based on binder value of two molecules
	character(len=256) ::   HBONDFileName
	integer*8, intent(IN) ::  binderi,&
							binderj
	
	if(binderi == binderj) then
		if (binderi > 3) then
			write(HBONDFileName, "(A,I1,A)") "HBONDS_dL", binderi-4, "_intra.bonds"
		else
			write(HBONDFileName, "(A,I1,A)") "HBONDS_uL", binderi, "_intra.bonds"
		end if
	else 
		if (binderi > 3) then
			write(HBONDFileName, "(A,I1,A)") "HBONDS_dL", binderi-4, "_inter.bonds"
		else
			write(HBONDFileName, "(A,I1,A)") "HBONDS_uL", binderi, "_inter.bonds"
		end if
	end if
end function HBONDFileName

function HistogramFileName(binderi)
	!addition to deployHistogramFiles()
	!returns name of the Histogram file based on binder value of two molecules
	character(len=256) ::   HistogramFileName
	integer*8, intent(IN) ::  binderi
	
	if(binderi < 8) then
		write(HistogramFileName, "(A,I1,A)") "Histogram_dL", binderi, ".3dhist"
	else 
		write(HistogramFileName, "(A,I1,A)") "Histogram_uL", binderi-8, ".3dhist"
	end if
	
end function HistogramFileName

end program Binder_P
