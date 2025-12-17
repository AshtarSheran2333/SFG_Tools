module SFG_FRAMES
	use SFG_TYPES
	implicit none
	
	private
		integer ::                                          noh = 0
		
		logical ::                                          hydroxyls = .false.
		
		character(10) ::                                    XNAME = "AL"
		
	public                                                  frame_reader,&
															trr_frame_reader,&
															gro_frame_reader,&
															xyz_frame_reader,&
															set_hydroxyl_metal_name
   
	!abstract type frame_reader to operate with pos+vel input files
	type, abstract :: frame_reader
		
		type(water), dimension(:), allocatable ::           molecule
		
		integer, private ::                                 file1,&
															file2
		
	contains
		procedure(open_), deferred ::                       open_file
		procedure(read_), deferred ::                       read_frame
		procedure(skip_), deferred ::                       skip_frame
		procedure(rewind_), deferred ::                     rewind_file
		procedure(is_hydroxyl_), deferred ::                is_hydroxyl
		procedure(close_), deferred ::                      close_file
		procedure(is_open_), deferred ::                    is_open 
		
	end type frame_reader
	
	abstract interface
	
		subroutine read_(this)
			import frame_reader
			
			class(frame_reader) ::                          this
		end subroutine read_
		!______________________________________________________________________
		subroutine skip_(this)
			import frame_reader
			
			class(frame_reader) ::                          this
		end subroutine skip_
		!______________________________________________________________________
		subroutine rewind_(this)
			import frame_reader
			
			class(frame_reader) ::                          this
		end subroutine rewind_
		!______________________________________________________________________
		subroutine open_(this, filename, filename2)
			import frame_reader
			
			class(frame_reader) ::                          this
			
			character(*), intent(IN) ::                     filename
			
			character(*), optional, intent(IN) ::           filename2
		end subroutine open_
		!______________________________________________________________________        
		logical function is_hydroxyl_(this)
			import frame_reader
			
			class(frame_reader) ::                          this
		end function
		!______________________________________________________________________        
		subroutine close_(this)
			import frame_reader
			
			class(frame_reader) ::                          this
		end subroutine
		!______________________________________________________________________
		logical function is_open_(this)
			import frame_reader
			
			class(frame_reader) ::                          this
		end function is_open_

		
	end interface
	
	!extedned frame_reader for trr files
	type, extends(frame_reader) :: trr_frame_reader
		
		integer::                                           natoms
		
	contains
	
		procedure ::                                        open_file => trr_open_file
		procedure ::                                        read_frame => trr_read_frame
		procedure ::                                        skip_frame => trr_skip_frame
		procedure ::                                        rewind_file => trr_rewind_file
		procedure ::                                        is_hydroxyl => trr_is_hydroxyl
		procedure ::                                        close_file => trr_close_file
		procedure ::                                        is_open => trr_is_open
		
	end type trr_frame_reader

	!extended frame_reader for gro files
	type, extends(frame_reader) :: gro_frame_reader
		
	contains
	
		procedure ::                                        open_file => gro_open_file
		procedure ::                                        read_frame => gro_read_frame
		procedure ::                                        skip_frame => gro_read_frame
		procedure ::                                        rewind_file => gro_rewind_file
		procedure ::                                        is_hydroxyl => gro_is_hydroxyl
		procedure ::                                        close_file => gro_close_file
		procedure ::                                        is_open => gro_is_open
		
	end type gro_frame_reader
	
	!extedned frame_reader for xyz files
	type, extends(frame_reader) :: xyz_frame_reader
		
	contains
		procedure ::                                        open_file => xyz_open_file
		procedure ::                                        read_frame => xyz_read_frame
		procedure ::                                        skip_frame => xyz_skip_frame
		procedure ::                                        rewind_file => xyz_rewind_file
		procedure ::                                        is_hydroxyl => xyz_is_hydroxyl
		procedure ::                                        close_file => xyz_close_file
		procedure ::                                        is_open => xyz_is_open
		
	end type xyz_frame_reader
	
contains

		
	subroutine set_hydroxyl_metal_name(name)
		character(10), intent(IN) :: name
		XNAME = adjustl(name)
	end subroutine set_hydroxyl_metal_name
	
	
	!################################trr frame_reader procedures########################################    
	
	subroutine trr_open_file(this, filename, filename2)
	
		class(trr_frame_reader) ::                          this
		
		character(*), intent(IN) ::                         filename
		character(*), optional, intent(IN) ::               filename2
		
		integer ::                                          natoms,&
															sumcheck = 0,&
															i,&
															int_num,&
															ierr
		
		character(len = 255) ::                             firstLine
		
		character(len = 10) ::                              atomtype,&
															atomtype1
		
		print*, "Opening *.trr file..."
		
		!if memory is allocated deallocate it
		if (allocated(this%molecule)) then
			deallocate(this%molecule)
		end if
		
		noh = 0
		
		open(newunit = this%file1, file = filename, form = "unformatted", access = 'stream', convert = 'big_endian', status = 'old', iostat = ierr)
		if(ierr .ne. 0) then
			print*, "unable to open ", trim(filename), " file"
			stop
		end if
		open(newunit = this%file2, file = filename2, iostat = ierr)
		if(ierr .ne. 0) then
			print*, "unable to open ", trim(filename2), " file"
			stop
		end if
		
		!first get info if the file contains hydroxyls and how much of them are connected to a metal atom
		read(this%file2,'(A)', iostat = ierr) firstline
		if(ierr .ne. 0) then
			print*, "corrupted *.gro file"
			stop
		end if
		
		read(this%file2,*, iostat = ierr) natoms
		if(ierr .ne. 0) then
			print*, "corrupted *.gro file"
			stop
		end if
				
		read(this%file2,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype, atomtype
		if(ierr .ne. 0) then
			print*, "corrupted *.gro file"
			stop
		end if
		
		if(adjustl(atomtype) == XNAME) then !in this case we opened a hydroxyl file
			hydroxyls = .true.
		
			read(this%file2,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype, atomtype
			if(ierr .ne. 0) then
				print*, "corrupted *.gro file"
				stop
			end if
					
			read(this%file2,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype1, atomtype1
			if(ierr .ne. 0) then
				print*, "corrupted *.gro file"
				stop
			end if
			
			do while((atomtype .NE. atomtype1) .AND. ((adjustl(atomtype) .NE. XNAME) .AND. (adjustl(atomtype1) .NE. XNAME)))
				noh = noh + 1
				read(this%file2,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype, atomtype
				if(ierr .ne. 0) then
					print*, "corrupted *.gro file"
					stop
				end if
						
				read(this%file2,"(i5,2a5,i5,3f8.3,3f8.4)") int_num, atomtype1, atomtype1
				if(ierr .ne. 0) then
					print*, "corrupted *.gro file"
					stop
				end if
			end do
			
			!rewind file, reopen it and count metals so we know that we need to allocate nsi+nsi*noh*2 "molecules"
			!also we can check if the number corresponds to natoms at the start of the file
			
			rewind(this%file2)
			
			read(this%file2,'(A)', iostat = ierr) firstline
			if(ierr .ne. 0) then
				print*, "corrupted *.gro file"
				stop
			end if
			
			read(this%file2,*, iostat = ierr) natoms
			if(ierr .ne. 0) then
				print*, "corrupted *.gro file"
				stop
			end if
			
			do i=1,natoms
				read(this%file2,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype, atomtype
				if(ierr .ne. 0) then
					print*, "corrupted *.gro file"
					stop
				end if
				
				if(adjustl(atomtype) == XNAME) then
					sumcheck = sumcheck+1
				end if
			end do
			
			if(sumcheck + sumcheck*noh*2 .EQ. natoms) then
				print*, "OK - gro file contains only hydroxyls"
				print*, "number of metal atoms: ", sumcheck
				print*, "number of O-H groups: ", noh
				print*, ""
				allocate(this%molecule(sumcheck*noh))
			else
				print*, "Input file structure does not contain only hydroxyls"
				stop !stopping the program because it is meaningless with bad input
			end if  
			
		else !in this case we opened just water file -> count oxygens and allocate number of oxygens "molecules"
			if(adjustl(atomtype) .EQ. "OW") then
				sumcheck = sumcheck+1
			end if
			
			do i=2,natoms
				read(this%file2,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype, atomtype
				if(ierr .ne. 0) then
					print*, "corrupted *.gro file"
					stop
				end if
				
				if(adjustl(atomtype) .EQ. "OW") then
					sumcheck = sumcheck+1
				end if
			end do
			
			if(3*sumcheck .EQ. natoms) then
				print*, "OK - gro file contains only waters"
				print*, "number of water molecules: ", sumcheck
				print*, ""
				allocate(this%molecule(sumcheck))
			else
				print*, "Input file structure does not contain only water molecules"
				stop !stopping the program because it is meaningless with bad input
			end if  
		end if
		
		!at this point we have allocated molecules that we need, we can close the gro file and never open it again
		this%natoms = natoms
		close(this%file2)
		
	end subroutine trr_open_file
	
	!__________________________________________________________________________
	
	subroutine trr_read_frame(this)
	
		class(trr_frame_reader) ::                          this
		
		real(8), dimension(3) ::                            metal
		
		character, dimension(:), allocatable ::             cw
		
		real(4) ::                                          read_data
		
		integer(4) ::                                       skip,&
															natoms
		
		integer ::                                          i,&
															j,&
															k,&
															read_status,&
															index
		
		character, dimension(12) ::                         control = "GMX_trn_file"
		
		logical ::                                          forcesFlag = .false.
		
	!verifying filetype
		
		read(this%file1) skip
	
		if(skip /= 1993) then
			print*, "corrupted TRR file"
			stop
		end if
	
		read(this%file1) skip
		allocate(cw(skip-1))

	!skip 32 bits
		read(this%file1) skip
	
	!check if cw contains "GMX_trn_file"
		read(this%file1) cw

	!skipping 7x 32bits (previously skipping 10x 32bits 0,9)
		do i=0, 6
			read(this%file1) skip
		end do
	!reading positions flag if not 0 positions are present
		read(this%file1) skip
	!reading velocities flag if not 0 velocities are present
		read(this%file1) skip
	!reading forces flag if not 0 forces are present
		read(this%file1) skip
		if(skip .NE. 0) then
			forcesFlag = .true.
		end if
	!reading data
	!read natoms
		read(this%file1) natoms
	!read step number
		read(this%file1) skip
	!skip 32 bits
		read(this%file1) skip
	!read simulation time
		read(this%file1) skip
	!skip 32 bits
		read(this%file1) skip
	!read box dimension X
		read(this%file1) skip
	!skip 3x 32 bits
		read(this%file1) skip
		read(this%file1) skip
		read(this%file1) skip
	!read box dimension Y
		read(this%file1) skip
	!skip 3x 32 bits
		read(this%file1) skip
		read(this%file1) skip
		read(this%file1) skip
	!read box dimension Z
		read(this%file1) skip
		
	!reading positions data are multiplied by 10 -> conversion from nm to Angstroms
		
		if(hydroxyls) then
		! ... read metal-O-H as water molecule so metal is allways H1, in case of multihydroxyls remember metal and copy it...
			index = 1
			
			do i=1, this%natoms, (noh*2)+1
					
					!reads metal position
					do j = 1, 3  
						read(this%file1) read_data
						metal(j) = 10*read_data
					end do
					
					do k = 1, noh
					
						!assign metal to each hydroxyl
						this%molecule(index)%h2%position = metal
						
						!read O position
						do j = 1, 3
							read(this%file1) read_data
							this%molecule(index)%o%position(j) = 10*read_data
						end do
						!read H position
						do j = 1, 3
							read(this%file1) read_data
							this%molecule(index)%h1%position(j) = 10*read_data
						end do
						!move to next metal-O-H index
						index = index + 1
					end do
			end do
			
		else
		! ... read H-O-H one after each other...
			do i=1, this%natoms, 3

					do j = 1, 3  
						read(this%file1) read_data
						this%molecule(ceiling(real(i)/3))%o%position(j) = 10*read_data
					end do
			
					do j = 1, 3
						read(this%file1) read_data
						this%molecule(ceiling(real(i)/3))%h1%position(j) = 10*read_data
					end do
			
					do j = 1, 3
						read(this%file1) read_data
						this%molecule(ceiling(real(i)/3))%h2%position(j) = 10*read_data
					end do
			end do
		end if
	!reading velocities multiplying by factor = 0.01 -> conversion from nm/ps to A/fs 
		if(hydroxyls) then 
			index = 1
		
			do i=1, this%natoms, (noh*2)+1
					
				!reads metal velocity
				do j = 1, 3  
					read(this%file1) read_data
					metal(j) = 0.01*read_data
				end do
					
				do k = 1, noh
					
					!assign metal to each hydroxyl
					this%molecule(index)%h1%velocity = metal
						
					!read O velocity
					do j = 1, 3
						read(this%file1) read_data
						this%molecule(index)%o%velocity(j) = 0.01*read_data
					end do
					!read H velocity
					do j = 1, 3
						read(this%file1) read_data
						this%molecule(index)%h2%velocity(j) = 0.01*read_data
					end do
					!move to next metal-O-H index
					index = index + 1
				end do
			end do
		else
			do i=1, this%natoms, 3

				do j = 1, 3  
					read(this%file1) read_data
					this%molecule(ceiling(real(i)/3))%o%velocity(j) = 0.01*read_data
				end do
			
				do j = 1, 3
					read(this%file1) read_data
					this%molecule(ceiling(real(i)/3))%h1%velocity(j) = 0.01*read_data
				end do

				do j = 1, 3
					read(this%file1) read_data
					this%molecule(ceiling(real(i)/3))%h2%velocity(j) = 0.01*read_data
				end do
			end do
		end if
		
		if(forcesFlag .eqv. .true.) then !if forces skip forces...
			do i = 1, this%natoms
				do j = 1,3
					read(this%file1) read_data
				end do
			end do
		end if

	end subroutine trr_read_frame

	!__________________________________________________________________________
	
	subroutine trr_skip_frame(this)
	
		class(trr_frame_reader) ::                          this
		
		character, dimension(:), allocatable ::             cw
		
		integer(4) ::                                       skip  
		
		integer ::                                          i,&
															ierr,&
															natoms 
		
		character, dimension(12) ::                         control = "GMX_trn_file"
		
		logical ::                                          posFlag = .false.,&
															velFlag = .false.,&
															forcesFlag = .false.
		
		!verifying file type
		read(this%file1, IOSTAT = ierr) skip
		if(ierr .ne. 0) then
			print*, "corrupted *.trr file"
			stop
		end if
	
		if(skip /= 1993) then
			print*, "corrupted *.trr file"
			stop
		end if
	
		read(this%file1, iostat = ierr) skip
		if(ierr .ne. 0) then
			print*, "corrupted *.trr file"
			stop
		end if
		
		allocate(cw(skip-1))
	!skip 32 bits
		read(this%file1, iostat = ierr) skip
		if(ierr .ne. 0) then
			print*, "corrupted *.trr file"
			stop
		end if
		
	!check if cw contains "GMX_trn_file"
		read(this%file1, iostat = ierr) cw
		if(ierr .ne. 0) then
			print*, "corrupted *.trr file"
			stop
		end if
		
		if(ALL(cw /= control)) then
			print*, "corrupted TRR file"
			stop
		end if
		
	!skipping 7x 32bits
		do i=0, 6
			read(this%file1, iostat = ierr) skip
			if(ierr .ne. 0) then
				print*, "corrupted *.trr file"
				stop
			end if
		end do
	!data flags
	!pos
		read(this%file1, iostat = ierr) skip
		if(ierr .ne. 0) then
			print*, "corrupted *.trr file"
			stop
		end if
		
		if(skip .ne. 0) then
			posFlag = .true.
		end if
	!vel
		read(this%file1, iostat = ierr) skip
		if(ierr .ne. 0) then
			print*, "corrupted *.trr file"
			stop
		end if
		
		if(skip .ne. 0) then
			velFlag = .true.
		end if
	!forces
		read(this%file1, iostat = ierr) skip
		if(ierr .ne. 0) then
			print*, "corrupted *.trr file"
			stop
		end if
		
		if(skip .ne. 0) then
			forcesFlag = .true.
		end if
		
		if((posFlag .eqv. .false.) .or. (velFlag .eqv. .false.)) then
			print*, "*.trr File does not contain both positions and velocities"
		end if
		
	!read natoms
		read(this%file1, iostat = ierr) natoms
		if(ierr .ne. 0) then
			print*, "corrupted *.trr file"
			stop
		end if
		
		if(natoms /= this%natoms) then
			print*, "*.trr file number of atoms does not match gro file number of atoms"
			stop
		end if
	!skipping 13x 32bits
		do i=0, 12
			read(this%file1, iostat = ierr) skip
			if(ierr .ne. 0) then
				print*, "corrupted *.trr file"
				stop
			end if
		end do
	!end of header
	!skipping positions and velocities and forces
		
		if(posFlag .eqv. .true.) then
			do i=1, 3*natoms
				read(this%file1, iostat = ierr) skip
				if(ierr .ne. 0) then
					print*, "corrupted *.trr file"
					stop
				end if
			end do
		end if
		
		if(velFlag .eqv. .true.) then
			do i=1, 3*natoms
				read(this%file1, iostat = ierr) skip
				if(ierr .ne. 0) then
					print*, "corrupted *.trr file"
					stop
				end if
			end do
		end if
		
		if(forcesFlag .eqv. .true.) then
			do i=1, 3*natoms
				read(this%file1, iostat = ierr) skip
				if(ierr .ne. 0) then
					print*, "corrupted *.trr file"
					stop
				end if
			end do
		end if
		
	end subroutine trr_skip_frame
	
	!__________________________________________________________________________
	
	logical function trr_is_hydroxyl(this)
	
			!import frame_reader
			
			class(trr_frame_reader) ::                      this
			
			if(hydroxyls) then
				trr_is_hydroxyl = .true.
			else
				trr_is_hydroxyl = .false.
			end if
			
	end function trr_is_hydroxyl
	
	!__________________________________________________________________________
	
	logical function trr_is_open(this)
	
			!import frame_reader
			
			class(trr_frame_reader) ::                      this
			inquire( unit = this%file1, opened = trr_is_open )
			
	end function trr_is_open
	
	!__________________________________________________________________________
	
	subroutine trr_close_file(this)
	
			!import frame_reader
			
			class(trr_frame_reader) ::                      this
			
			close(this%file1)
			
	end subroutine trr_close_file
	
	
	!__________________________________________________________________________
	
	subroutine trr_rewind_file(this)
	
		!import frame_reader
			
		class(trr_frame_reader) ::                      this
			
		rewind(this%file1)
			
	end subroutine trr_rewind_file
	
	!################################gro frame_reader procedures########################################  
	
	subroutine gro_open_file(this, filename, filename2)
	
		class(gro_frame_reader) ::                          this
		
		character(*), intent(IN) ::                         filename
		
		character(*), optional, intent(IN) ::               filename2
		
		integer ::                                          natoms,&
															sumcheck = 0,&
															i,&
															int_num,&
															ierr
		
		character(len = 255) ::                             firstLine
		
		character(len = 10) ::                              atomtype,&
															atomtype1
		
		print*, "Opening *.gro file..."
		
		
		!if something is allocated, deallocate 
		if (allocated(this%molecule)) then
			deallocate(this%molecule)
		end if
	
		noh = 0
		
		open(newunit = this%file1, file = filename, iostat = ierr)
		if(ierr .ne. 0) then
			print*, "unable to open ", trim(filename), " file"
			stop
		end if
		
		!first get info if the file contains hydroxyls and how much of them are connected to a metal atom
		read(this%file1,'(A)', iostat = ierr) firstline
		if(ierr .ne. 0) then
			print*, "corrupted *.gro file"
			stop
		end if
		
		read(this%file1,*, iostat = ierr) natoms
		if(ierr .ne. 0) then
			print*, "corrupted *.gro file"
			stop
		end if
		
		read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype, atomtype
		if(ierr .ne. 0) then
			print*, "corrupted *.gro file"
			stop
		end if
		
		if(adjustl(atomtype) == XNAME) then !in this case we opened a hydroxyl file
			hydroxyls = .true.
		
			read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype, atomtype
			if(ierr .ne. 0) then
				print*, "corrupted *.gro file"
				stop
			end if
			
			read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype1, atomtype1
			if(ierr .ne. 0) then
				print*, "corrupted *.gro file"
				stop
			end if
			
			do while((atomtype .NE. atomtype1) .AND. ((adjustl(atomtype) .NE. XNAME) .AND. (adjustl(atomtype1) .NE. XNAME)))
				noh = noh + 1
				read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype, atomtype
				if(ierr .ne. 0) then
					print*, "corrupted *.gro file"
					stop
				end if
				
				read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype1, atomtype1
				if(ierr .ne. 0) then
					print*, "corrupted *.gro file"
					stop
				end if
			end do
			
			!close file, reopen it and count metals so we know that we need to allocate nsi*noh "molecules"
			!also we can check if the number corresponds to natoms at the start of the file
			
			rewind(this%file1)
			
			read(this%file1,'(A)', iostat = ierr) firstline
			if(ierr .ne. 0) then
				print*, "corrupted *.gro file"
				stop
			end if
			
			read(this%file1,*, iostat = ierr) natoms
			if(ierr .ne. 0) then
				print*, "corrupted *.gro file"
				stop
			end if
			
			do i=1,natoms
				read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype, atomtype
				if(ierr .ne. 0) then
					print*, "corrupted *.gro file"
					stop
				end if
				
				if(adjustl(atomtype) == XNAME) then 
					sumcheck = sumcheck+1 
				end if
			end do
			
			if(sumcheck + sumcheck*noh*2 .EQ. natoms) then
				print*, "OK - gro file contains only hydroxyls"
				print*, "number of metal atoms: ", sumcheck
				print*, "number of O-H groups: ", noh
				print*, ""
				allocate(this%molecule(sumcheck*noh))
			else
				print*, "Input file does not contain only hydroxyls"
				stop !stopping the program because it is meaningless with bad input
			end if  
			
		else !in this case we opened just water file -> count oxygens and allocate number of oxygens "molecules"
			if(adjustl(atomtype) .EQ. "OW") then
				sumcheck = sumcheck+1
			end if
			
			do i=2,natoms
				read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) int_num, atomtype, atomtype
				if(ierr .ne. 0) then
					print*, "corrupted *.gro file"
					stop
				end if
				if(adjustl(atomtype) .EQ. "OW") then
					sumcheck = sumcheck+1
				end if
			end do
			
			if(3*sumcheck .EQ. natoms) then
				print*, "OK - gro file contains only waters"
				print*, "number of water molecules: ", sumcheck
				print*, ""
				allocate(this%molecule(sumcheck))
			else
				print*, "Input file does not contain only water molecules"
				stop !stopping the program because it is meaningless with bad input
			end if  
		end if
		
		!at this point we have allocated molecules that we need, reopen gro file to start reading from beginning
		rewind(this%file1)
		
	end subroutine gro_open_file
	
	!__________________________________________________________________________
	
	subroutine gro_read_frame(this)
		class(gro_frame_reader) ::                          this
		
		type(atom) ::                                       metal
		
		real(8) ::                                          x,&
															y,&
															z,&
															vx,&
															vy,&
															vz
		
		integer ::                                          i,&
															j,&
															atomNumber,&
															natoms,&
															index,&
															ierr
		
		character(len = 255) ::                             firstLine
		
		character(len = 5) ::                               atomType
		
		character(len = 5) ::                               atomName
		
		logical ::                                          velocities=.true.

		!look for "t="
		!reads whole first line as a string
		read(this%file1,'(A)') firstLine
		!reading number of atoms in this frame
		read(this%file1,*) natoms
		
		!sorting the gro frame
		
		if(hydroxyls) then
			index = 1
			do i = 1, natoms, 2*noh+1
				
				
				if(i == 1) then !check if velocities are present in the first line
					read(this%file1,'(A)') firstline
					read(firstline, "(i5,2a5,i5,3f8.3,3f8.4)", IOSTAT = ierr) atomNumber, atomType, atomName, atomNumber, x, y, z, vx, vy, vz
					if(ierr .ne. 0) then
						read(firstline, "(i5,2a5,i5,3f8.3,3f8.4)", IOSTAT = ierr) atomNumber, atomType, atomName, atomNumber, x, y, z
						velocities = .false.
					end if
					
				else 
					if(velocities .eqv. .true.) then
						read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)") atomNumber, atomType, atomName, atomNumber, x, y, z, vx, vy, vz
					else
						read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)") atomNumber, atomType, atomName, atomNumber, x, y, z
					end if
				end if
				
					!save metal positions from nm to A, velocities from nm/ps to A/fs
					metal%position(1) = 10*x
					metal%position(2) = 10*y
					metal%position(3) = 10*z
					metal%velocity(1) = 0.01*vx
					metal%velocity(2) = 0.01*vy
					metal%velocity(3) = 0.01*vz
				
				do j = 1, noh
					!assign metal
					this%molecule(index)%h2 = metal
					
					!read O
					if(velocities .eqv. .true.) then
						read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)") atomNumber, atomType, atomName, atomNumber, x, y, z, vx, vy, vz
					else
						read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)") atomNumber, atomType, atomName, atomNumber, x, y, z
					end if
					this%molecule(index)%o%position(1) = 10*x
					this%molecule(index)%o%position(2) = 10*y
					this%molecule(index)%o%position(3) = 10*z
					this%molecule(index)%o%velocity(1) = 0.01*vx
					this%molecule(index)%o%velocity(2) = 0.01*vy
					this%molecule(index)%o%velocity(3) = 0.01*vz    
					
					!read H
					if(velocities .eqv. .true.) then
						read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)") atomNumber, atomType, atomName, atomNumber, x, y, z, vx, vy, vz
					else
						read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)") atomNumber, atomType, atomName, atomNumber, x, y, z
					end if
					this%molecule(index)%h1%position(1) = 10*x
					this%molecule(index)%h1%position(2) = 10*y
					this%molecule(index)%h1%position(3) = 10*z
					this%molecule(index)%h1%velocity(1) = 0.01*vx
					this%molecule(index)%h1%velocity(2) = 0.01*vy
					this%molecule(index)%h1%velocity(3) = 0.01*vz
					
					index = index + 1
				end do
			end do
		else
			index = 0
			
			do i = 1, natoms
				
				if(i == 1) then !check if velicities are present in the first line
					read(this%file1,'(A)') firstline
					read(firstline, "(i5,2a5,i5,3f8.3,3f8.4)", IOSTAT = ierr) atomNumber, atomType, atomName, atomNumber, x, y, z, vx, vy, vz
					if(ierr .ne. 0) then
						read(firstline, "(i5,2a5,i5,3f8.3,3f8.4)", IOSTAT = ierr) atomNumber, atomType, atomName, atomNumber, x, y, z
						velocities = .false.
					end if
					
				else 
					if(velocities .eqv. .true.) then
						read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)") atomNumber, atomType, atomName, atomNumber, x, y, z, vx, vy, vz
					else
						read(this%file1,"(i5,2a5,i5,3f8.3,3f8.4)") atomNumber, atomType, atomName, atomNumber, x, y, z
					end if
				end if
		
				select case(adjustl(atomName))
		
				case("OW")
					index = index + 1
					this%molecule(index)%o%position(1) = 10*x
					this%molecule(index)%o%position(2) = 10*y
					this%molecule(index)%o%position(3) = 10*z
					this%molecule(index)%o%velocity(1) = 0.01*vx
					this%molecule(index)%o%velocity(2) = 0.01*vy
					this%molecule(index)%o%velocity(3) = 0.01*vz
				case("HW1")
					this%molecule(index)%h1%position(1) = 10*x
					this%molecule(index)%h1%position(2) = 10*y
					this%molecule(index)%h1%position(3) = 10*z
					this%molecule(index)%h1%velocity(1) = 0.01*vx
					this%molecule(index)%h1%velocity(2) = 0.01*vy
					this%molecule(index)%h1%velocity(3) = 0.01*vz           
				case("HW2")
					this%molecule(index)%h2%position(1) = 10*x
					this%molecule(index)%h2%position(2) = 10*y
					this%molecule(index)%h2%position(3) = 10*z
					this%molecule(index)%h2%velocity(1) = 0.01*vx
					this%molecule(index)%h2%velocity(2) = 0.01*vy
					this%molecule(index)%h2%velocity(3) = 0.01*vz           
				case default
			
				end select
			end do
		end if
			
		!reading box dimensions
		read(this%file1,*) x,y,z
	
	end subroutine gro_read_frame
	
	!__________________________________________________________________________
	
	subroutine gro_skip_frame(this)
	
		class(gro_frame_reader) ::                          this
		
		integer ::                                          i,&
															NTOT,&
															atomNumber,&
															ierr
		
		character(len = 255) ::                             firstLine
		
		character(len = 5) ::                               atomType
		
		character(len = 5) ::                               atomName
		
		real(kind=8) ::                                           x,&
															y,&
															z
		
		!reads whole first line as a string
		read(this%file1,'(A)', iostat = ierr) firstLine
		if(ierr .ne. 0) then
			print*, "ERROR: *.gro file is corrupted"
			stop
		end if
		
		if(index(firstLine, "t=") .EQ. 0) then
			print*, "ERROR: *.gro file is corrupted"
		end if
		
		read(this%file1,*, iostat = ierr) NTOT
		if(ierr .ne. 0) then
			print*, "ERROR: *.gro file is corrupted"
			stop
		end if
		
		do i = 0, NTOT
			read(this%file1, "(i5,2a5,i5,3f8.3,3f8.4)", iostat = ierr) atomNumber, atomType, atomName, atomNumber, x, y, z
			
			if(ierr .ne. 0) then
				print*, "ERROR: *.gro file is corrupted"
				stop
			end if
			
		end do
		
	end subroutine gro_skip_frame
		
	!__________________________________________________________________________
	
	logical function gro_is_hydroxyl(this)
	
		!import frame_reader
			
		class(gro_frame_reader) ::                      this
			
		if(hydroxyls) then
			gro_is_hydroxyl = .true.
		else
			gro_is_hydroxyl = .false.
		end if
			
	end function gro_is_hydroxyl
		
	!__________________________________________________________________________
	
	logical function gro_is_open(this)
	
			!import frame_reader
			
			class(gro_frame_reader) ::                      this
			inquire( unit = this%file1, opened = gro_is_open )
			
	end function gro_is_open
	
	!__________________________________________________________________________
	
	subroutine gro_close_file(this)
	
		!import frame_reader
	
		class(gro_frame_reader) ::                      this
			
		close(this%file1)
			
	end subroutine gro_close_file
	
	 !__________________________________________________________________________
	
	subroutine gro_rewind_file(this)
	
		!import frame_reader
			
		class(gro_frame_reader) ::                      this
			
		rewind(this%file1)
			
	end subroutine gro_rewind_file
	
	!################################xyz frame_reader procedures######################################## 
	
	subroutine xyz_open_file(this, filename, filename2)
	
		class(xyz_frame_reader) ::                          this
		
		character(*), intent(IN) ::                         filename
		
		character(*), optional, intent(IN) ::               filename2
		
		integer ::                                          natoms,&
															sumcheck = 0,&
															i,&
															ierr
		
		character(len = 255) ::                             firstLine
		
		character(len = 10) ::                              atomtype, atomtype1
		
		print*, "Opening *.xyz files..."
		
		!opens both xyz files
		open (newunit = this%file1, file = filename, status = 'old', iostat = ierr)
		if(ierr .ne. 0) then
			print*, "unable to open ", trim(filename), " file"
			stop
		end if
		
		open (newunit = this%file2, file = filename2, status = 'old', iostat = ierr)
		if(ierr .ne. 0) then
			print*, "unable to open ", trim(filename2), " file"
			stop
		end if
		
		if (allocated(this%molecule)) then
			deallocate(this%molecule)
		end if
	
		noh = 0
		
		!first get knowledge about if the file contains hydroxyls and how much of them are connected to a metal atom
		read(this%file1,*, iostat = ierr) natoms
		if(ierr .ne. 0) then
			print*, "corrupted *.xyz file"
			stop
		end if
		
		read(this%file1,'(A)', iostat = ierr) firstline
		if(ierr .ne. 0) then
			print*, "corrupted *.xyz file"
			stop
		end if
		
		read(this%file1,*, iostat = ierr) atomtype
		if(ierr .ne. 0) then
			print*, "corrupted *.xyz file"
			stop
		end if
		
		if(atomtype == XNAME) then !in this case we opened a hydroxyl file
			hydroxyls = .true.
		
			read(this%file1,*, iostat = ierr) atomtype
			if(ierr .ne. 0) then
				print*, "corrupted *.xyz file"
				stop
			end if
			read(this%file1,*, iostat = ierr) atomtype1
			if(ierr .ne. 0) then
				print*, "corrupted *.xyz file"
				stop
			end if
			
			do while((atomtype .NE. atomtype1) .AND. ((atomtype .NE. XNAME) .AND. (atomtype1 .NE. XNAME)))
				noh = noh + 1
				read(this%file1,*, iostat = ierr) atomtype
				if(ierr .ne. 0) then
					print*, "corrupted *.xyz file"
					stop
				end if
				
				read(this%file1,*, iostat = ierr) atomtype1
				if(ierr .ne. 0) then
					print*, "corrupted *.xyz file"
					stop
				end if
			end do
			
			!close file, reopen it and count metals so we know that we need to allocate nmetal*noh "molecules"
			!also we can check if the number corresponds to natoms at the start of the file
			
			rewind(this%file1)
			
			read(this%file1,*, iostat = ierr) natoms
			if(ierr .ne. 0) then
				print*, "corrupted *.xyz file"
				stop
			end if
			
			read(this%file1,'(A)', iostat = ierr) firstline
			if(ierr .ne. 0) then
				print*, "corrupted *.xyz file"
				stop
			end if
			
			do i=1,natoms
				read(this%file1,*, iostat = ierr) atomtype
				if(ierr .ne. 0) then
					print*, "corrupted *.xyz file"
					stop
				end if
				
				if(atomtype == XNAME) then 
					sumcheck = sumcheck+1 
				end if
			end do
			
			if(sumcheck + sumcheck*noh*2 .EQ. natoms) then
				print*, "OK - xyz file contains only hydroxyls"
				print*, "number of metal atoms: ", sumcheck
				print*, "number of O-H groups: ", noh
				print*, ""
				allocate(this%molecule(sumcheck*noh))
			else
				print*, "Input file does not contain only hydroxyls"
				stop !stopping the program because it is meaningless with bad input
			end if  
			
		else !in this case we opened just water file -> count oxygens and allocate number of oxygens "molecules"
			if(atomtype .EQ. "O") then
				sumcheck = sumcheck+1
			end if
			
			do i=2,natoms
				read(this%file1,*, iostat = ierr) atomtype
				if(ierr .ne. 0) then
					print*, "corrupted *.xyz file"
					stop
				end if
				
				if(atomtype .EQ. "O") then
					sumcheck = sumcheck+1
				end if
			end do
			
			if(3*sumcheck .EQ. natoms) then
				print*, "OK - xyz file contains only waters"
				print*, "number of water molecules: ", sumcheck
				print*, ""
				allocate(this%molecule(sumcheck))
			else
				print*, "Input file does not contain only water molecules"
				stop !stopping the program because it is meaningless with bad input
			end if  
		end if
		
		!at this point we have allocated molecules that we need, reopen xyz_position file to start reading from beginning
		rewind(this%file1)
		
	end subroutine xyz_open_file
	
	!__________________________________________________________________________
	
	subroutine xyz_read_frame(this)
	
		class(xyz_frame_reader) ::                          this
		
		type(atom) ::                                       metal
		
		real(8) ::                                          x,&
															y,&
															z,&
															vx,&
															vy,&
															vz
		
		integer ::                                          i,&
															j,&
															POS_NATOM,&
															VEL_NATOM,&
															index,&
															NH
		
		character(len = 3) ::                               pos_atomName,&
															vel_atomName
		
		!reading number of atoms in both files
		read(this%file1,*) POS_NATOM
		read(this%file2,*) VEL_NATOM

		!reads comment line in both xyz files
		read(this%file1,*)
		read(this%file2,*)
		
		!if input files contains hydroxyls !velocities are converted from Hartree atomic units to A/fs, positions are already in A...
		if(hydroxyls) then
			index = 1
			do i = 1, POS_NATOM, 2*noh+1
				read(this%file1,*) pos_atomName, x, y, z
				read(this%file2,*) vel_atomName, vx, vy, vz

					!save metal
					metal%position(1) = x
					metal%position(2) = y
					metal%position(3) = z
					metal%velocity(1) = vx*21.8769125400
					metal%velocity(2) = vy*21.8769125400
					metal%velocity(3) = vz*21.8769125400
				
				do j = 1, noh
					!assign metal
					this%molecule(index)%h2 = metal
					
					!read O
					read(this%file1,*) pos_atomName, x, y, z
					read(this%file2,*) vel_atomName, vx, vy, vz

					this%molecule(index)%o%position(1) = x
					this%molecule(index)%o%position(2) = y
					this%molecule(index)%o%position(3) = z
					this%molecule(index)%o%velocity(1) = vx*21.8769125400
					this%molecule(index)%o%velocity(2) = vy*21.8769125400
					this%molecule(index)%o%velocity(3) = vz*21.8769125400
					
					!read H
					read(this%file1,*) pos_atomName, x, y, z
					read(this%file2,*) vel_atomName, vx, vy, vz

					this%molecule(index)%h1%position(1) = x
					this%molecule(index)%h1%position(2) = y
					this%molecule(index)%h1%position(3) = z
					this%molecule(index)%h1%velocity(1) = vx*21.8769125400
					this%molecule(index)%h1%velocity(2) = vy*21.8769125400
					this%molecule(index)%h1%velocity(3) = vz*21.8769125400
					
					index = index + 1
				end do
			end do
		else
			index = 0
			NH = 0
			!reads both xyz files and sorts positions and velocities
			do i=1, POS_NATOM
				read(this%file1,*) pos_atomName, x, y, z
				read(this%file2,*) vel_atomName, vx, vy, vz
			
				select case(pos_atomName)
		
				case("O")
					index = index+1
					this%molecule(index)%o%position(1) = x
					this%molecule(index)%o%position(2) = y
					this%molecule(index)%o%position(3) = z
					this%molecule(index)%o%velocity(1) = vx*21.8769125400
					this%molecule(index)%o%velocity(2) = vy*21.8769125400
					this%molecule(index)%o%velocity(3) = vz*21.8769125400
				case("H")
					NH = NH+1
					if(mod(NH,2) == 1) then
						this%molecule(index)%h1%position(1) = x
						this%molecule(index)%h1%position(2) = y
						this%molecule(index)%h1%position(3) = z
						this%molecule(index)%h1%velocity(1) = vx*21.8769125400
						this%molecule(index)%h1%velocity(2) = vy*21.8769125400
						this%molecule(index)%h1%velocity(3) = vz*21.8769125400
					else
						this%molecule(index)%h2%position(1) = x
						this%molecule(index)%h2%position(2) = y
						this%molecule(index)%h2%position(3) = z
						this%molecule(index)%h2%velocity(1) = vx*21.8769125400
						this%molecule(index)%h2%velocity(2) = vy*21.8769125400
						this%molecule(index)%h2%velocity(3) = vz*21.8769125400
					endif
				case default
			
				end select
			end do
		end if
		
	end subroutine xyz_read_frame
	
	!__________________________________________________________________________
	
	subroutine xyz_skip_frame(this)
	
		class(xyz_frame_reader) ::                          this
		
		integer ::                                          i,&
															POS_NATOM,&
															VEL_NATOM,&
															ierr
		
		real(kind=8) ::                                           dummy_real
		character(len=3) ::                                 vel_dummy_string,&
															pos_dummy_string
		
		!reading number of atoms in both files
		read(this%file1,*, iostat = ierr) POS_NATOM
		if(ierr .ne. 0) then
			print*, "ERROR: position *.xyz file is corrupted"
			stop
		end if
		
		read(this%file2,*, iostat = ierr) VEL_NATOM
		if(ierr .ne. 0) then
			print*, "ERROR: velocity *.xyz file is corrupted"
			stop
		end if
		
		if(POS_NATOM /= VEL_NATOM) then
			   print*, ".xyz files do not match in number of atoms"
			   stop
		end if
		
		!reads comment line in both xyz files
		read(this%file1,*, iostat = ierr)
		if(ierr .ne. 0) then
			print*, "ERROR: position *.xyz file is corrupted"
			stop
		end if
		
		read(this%file2,*, iostat = ierr)
		if(ierr .ne. 0) then
			print*, "ERROR: velocity *.xyz file is corrupted"
			stop
		end if
		!reads both xyz files
		do i=1, POS_NATOM
			read(this%file1,*, iostat = ierr) pos_dummy_string, dummy_real, dummy_real, dummy_real
			if(ierr .ne. 0) then
				print*, "ERROR: position *.xyz file is corrupted"
				stop
			end if
			
			read(this%file2,*, iostat = ierr) vel_dummy_string, dummy_real, dummy_real, dummy_real
			if(ierr .ne. 0) then
				print*, "ERROR: velocity *.xyz file is corrupted"
				stop
			end if
			
			if(pos_dummy_string /= vel_dummy_string) then
				print*, "ERROR: *.xyz files does not match in atom structure"
				stop
			end if
		end do
		
	end subroutine xyz_skip_frame
	
	!__________________________________________________________________________
	
	logical function xyz_is_hydroxyl(this)
	
			!import frame_reader
			
			class(xyz_frame_reader) ::                      this
			
			if(hydroxyls) then
				xyz_is_hydroxyl = .true.
			else
				xyz_is_hydroxyl = .false.
			end if
			
	end function xyz_is_hydroxyl
		
	!__________________________________________________________________________
	
	logical function xyz_is_open(this)
	
			!import frame_reader
			
			class(xyz_frame_reader) ::                      this
			inquire( unit = this%file1, opened = xyz_is_open )
			inquire( unit = this%file2, opened = xyz_is_open )
			
	end function xyz_is_open
	
	!__________________________________________________________________________
	
	subroutine xyz_close_file(this)
	
		!import frame_reader
			
		class(xyz_frame_reader) ::                      this
			
		close(this%file1)
		close(this%file2)
			
	end subroutine xyz_close_file
	
	 !__________________________________________________________________________
	
	subroutine xyz_rewind_file(this)
	
		!import frame_reader
			
		class(xyz_frame_reader) ::                      this
			
		rewind(this%file1)
		rewind(this%file2)
			
	end subroutine xyz_rewind_file
	 
end module SFG_FRAMES
