module SFG_CORRELATION
	use, intrinsic :: iso_fortran_env
	use BINDER_FILE, only: group_binder_type
	use BOXDATA, only: boxdata_type
    use FRAME_READERS, only: current_frame_type
    use SFG_STRUCTURE, only: sfg_structure_group_type
    use DXDRZ_DB, only: dXdrz_db_type
	implicit none
    
	type correlation_function_type
        logical, dimension(:,:), allocatable :: binder_history !(SFG_UNITS, corrlen)
		real(real64), dimension(:,:), allocatable :: A_history, M_history !(SFG_UNITS, corrlen)
        real(real64), dimension(:), allocatable :: avA, avM !(SFG_UNITS)
        real(real64), dimension(:), allocatable :: correlation_function !(corrlen)
        integer(int64), dimension(:), allocatable :: norm !(corrlen)

        integer(int64) ::	corrlen,& !number of samples based on boxdata $corrlen
							t,& !time counter - used for the circular buffer
                            samples
        
        
    contains
		procedure, public :: init
        procedure, public :: calculate_step
        !procedure, public :: skip_step
    end type
    
    contains
    
    subroutine init(this, struct_group, boxdata)
        class(correlation_function_type), intent(inout) :: this
        type(sfg_structure_group_type), intent(in) :: struct_group
        type(boxdata_type), intent(in) :: boxdata

        !TODO could verify inputs e.g. struct_group has > 0 elements, ...

        !get rid of the old garbage if present
        if(allocated(this%binder_history)) deallocate(this%binder_history)
        if(allocated(this%A_history)) deallocate(this%A_history)
        if(allocated(this%M_history)) deallocate(this%M_history)
        if(allocated(this%avA)) deallocate(this%avA)
        if(allocated(this%avM)) deallocate(this%avM)
        if(allocated(this%correlation_function)) deallocate(this%correlation_function)
        if(allocated(this%norm)) deallocate(this%norm)
		
        !boxdata%CORRLEN (ps)
        !boxdata%DT (fs)
        this%corrlen = INT((1000 * boxdata%CORRLEN) / boxdata%DT) + 1
        this%t = 0 !nothing registered
        this%samples = 0 !used to get the average
        
        !allocate everything...
        allocate(this%binder_history(struct_group%n_elements, this%corrlen))
        allocate(this%A_history(struct_group%n_elements, this%corrlen))
        allocate(this%M_history(struct_group%n_elements, this%corrlen))
        allocate(this%avA(struct_group%n_elements))
        allocate(this%avM(struct_group%n_elements))
        allocate(this%correlation_function(this%corrlen))
        allocate(this%norm(this%corrlen))
        
        this%A_history = 0
        this%M_history = 0
        this%avA = 0
        this%avM = 0
        this%correlation_function = 0
        this%norm = 0
    end subroutine init

    subroutine calculate_step(this, current_frame, struct_group, binder, dXdrz_db, layer_selection, boxdata)
        class(correlation_function_type), intent(inout) :: this
        type(current_frame_type), intent(in) :: current_frame
        type(sfg_structure_group_type), intent(in) :: struct_group
        type(group_binder_type), intent(in) :: binder
        type(dXdrz_db_type), intent(in) :: dXdrz_db
        type(boxdata_type), intent(in) :: boxdata
        integer(int8), dimension(:), allocatable :: layer_selection

        integer(int64) :: nt1, nt0, timelag, m, l
        real(real64) :: weight
        real(real64), dimension(2) :: AM
        
        nt1 = mod(this%t, this%corrlen) + 1 !ringbuffer index
        this%t = this%t + 1
        !clear binder
        this%binder_history(:,nt1) = .false.
        
        !TODO make a function, OMP inside
        !go over all SFG units
        do m = 1, struct_group%n_elements
            !get am of m
            AM = struct_group%sfg_units(m)%get_AM(current_frame, dXdrz_db, boxdata)
            !append a m
            this%A_history(m, nt1) = AM(1)
            this%M_history(m, nt1) = AM(2)
            !update averages
            this%avA(m) = this%avA(m) + (AM(1) - this%avA(m))/this%t 
            this%avM(m) = this%avM(m) + (AM(2) - this%avM(m))/this%t 
            !subtract the best averages currently available
            this%A_history(m, nt1) = this%A_history(m,nt1) - this%avA(m)
            this%M_history(m, nt1) = this%M_history(m,nt1) - this%avM(m)

            !fill binder
            do l = 1, size(layer_selection)
                if(binder%layers(m) == layer_selection(l)) then
                    this%binder_history(m,nt1) = .true.
                    exit !this m was selected
                end if
            end do
        end do

        !TODO OMP
      	do timelag = 1, min(this%t,this%corrlen) !ramping the iterations from beginning...
            !the timelag actually goes from 0 to max_lag (corrlen - 1) indexing issues...
			this%norm(timelag) = this%norm(timelag) + 1
            !NT0 based on timelag, scan the whole history (wraparound of the ringbuffer)
            nt0 = mod(this%t - timelag, this%corrlen) + 1
            
			do m = 1, struct_group%n_elements
				! if the molecule is not in selected layer continue
				if(.not. (this%binder_history(m,nt0) .or. this%binder_history(m, nt1))) cycle
				
				if(this%binder_history(m,nt0) .and. this%binder_history(m,nt1)) then
					weight = 1.0_real64
				else
				! molecule M is present only in time nt1 or nt2
				! half weight for the term
                    weight = 0.5_real64
				end if
                
				this%correlation_function(timelag) = this%correlation_function(timelag) &
                        + (this%A_history(m, nt1)) * (this%M_history(m, nt0))
			end do
		end do
    
    end subroutine calculate_step

!    subroutine skip_step(this)
!        !fill binder
!    
!        !fill AM
!    end subroutine skip_step
    
	!binder in time(corrlen)
    !A, M (SFG_UNITS, corrlen)
    !avA, avM(SFG_UNITS)
    
	!if(mod(t,bd%self_skip) == 0) then 
	!	!$OMP PARALLEL DEFAULT(SHARED)
	!	!$OMP DO PRIVATE(timelag, nt0, m)
	!	! self correlation function
	!	do timelag = 0, min(t,bd%get_maxlag())
	!		selfnorm(timelag) = selfnorm(timelag) + 1
	!		nt0 = mod(t-timelag,corrlen)
	!		do m = 1, bd%NO
	!			! if the molecule is not in selected layer continue
	!			if(.not. (binder_in_time(m,nt0) .or. binder_in_time(m, nt1))) cycle
	!			
	!			! molecule M is present in both times nt1 and nt0...
	!			! full weight of the term
	!			if(binder_in_time(m,nt0) .and. binder_in_time(m,nt1)) then
	!				self_corr(timelag) = self_corr(timelag) + (Axx(m, nt1) - avAxx(m)) * (Mz(m, nt0)-avMz(m))
	!			else
	!			! molecule M is present only in time nt1 or nt2
	!			! half weight for the term
	!				self_corr(timelag) = self_corr(timelag) + 0.5d0 * (Axx(m, nt1) - avAxx(m)) * (Mz(m, nt0)-avMz(m))
	!			end if
	!		end do
	!	end do
	!	!$OMP END DO
	!	!$OMP END PARALLEL
	!end if
        
end module SFG_CORRELATION
