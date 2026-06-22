module SFG_CORRELATION
use, intrinsic :: iso_fortran_env
use BINDER_FILE, only: group_binder_type
use BOXDATA, only: boxdata_type
use FRAME_READERS, only: current_frame_type
use SFG_STRUCTURE, only: sfg_structure_group_type
use DXDRZ_DB, only: dXdrz_db_type
use SFG_UTILS, only: e_to_c, debye_to_ea
implicit none
    
type correlation_function_type
        logical, dimension(:,:), allocatable :: binder_history !(SFG_UNITS, corrlen)
        real(real64), dimension(:,:), allocatable :: A_history, M_history !(SFG_UNITS, corrlen)
        real(real64), dimension(:), allocatable :: avA, avM !(SFG_UNITS)
        real(real64), dimension(:), allocatable :: correlation_function !(corrlen)
        integer(int64), dimension(:), allocatable :: norm !(corrlen)

        integer(int64) ::corrlen,& !number of samples based on boxdata $corrlen
                            t, nt1 !time counter - used for the circular buffer
        
    contains
        procedure, public :: init
        procedure, public :: calculate_step
        procedure, public :: skip_step
        procedure, public :: get_normalized
        procedure, private :: fill_history
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
        integer(int8), dimension(:), allocatable :: layer_selection
        type(boxdata_type), intent(in) :: boxdata

        integer(int64) :: nt0, timelag, m
        real(real64) :: weight
        
        call this%fill_history(current_frame, struct_group, binder, dXdrz_db, layer_selection, boxdata)

        if(this%t == 1) return !in the first step, we would compute correlation of only zeroes, skip it

        !in this OMP loop, we are only writing to array(timelag), it is safe to share the whole instance of this
        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(this, struct_group) &
        !$OMP PRIVATE(timelag, nt0, m, weight)
        do timelag = 1, min(this%t,this%corrlen) !ramping the iterations from beginning...
            !the timelag actually goes from 0 to max_lag (corrlen - 1) indexing issues...
            this%norm(timelag) = this%norm(timelag) + 1
            !NT0 based on timelag, scan the whole history (wraparound of the ringbuffer)
            nt0 = mod(this%t - timelag, this%corrlen) + 1
            
            do m = 1, struct_group%n_elements
                ! if the molecule is not in selected layer continue
                if(.not. (this%binder_history(m,nt0) .or. this%binder_history(m, this%nt1))) cycle
            
                if(this%binder_history(m,nt0) .and. this%binder_history(m,this%nt1)) then
                    weight = 1.0_real64
                else
                    ! molecule M is present only in time nt1 or nt2
                    ! half weight for the term
                    weight = 0.5_real64
                end if
                
                this%correlation_function(timelag) = this%correlation_function(timelag) &
                        + (this%A_history(m, this%nt1)) * (this%M_history(m, nt0))
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine calculate_step

    subroutine fill_history(this, current_frame, struct_group, binder, dXdrz_db, layer_selection, boxdata)
        class(correlation_function_type), intent(inout) :: this
        type(current_frame_type), intent(in) :: current_frame
        type(sfg_structure_group_type), intent(in) :: struct_group
        type(group_binder_type), intent(in) :: binder
        type(dXdrz_db_type), intent(in) :: dXdrz_db
        integer(int8), dimension(:), allocatable :: layer_selection
        type(boxdata_type), intent(in) :: boxdata

        integer(int64) :: m, l
        real(real64), dimension(2) :: AM

        this%nt1 = mod(this%t, this%corrlen) + 1 !ringbuffer index
        this%t = this%t + 1
        !clear binder
        this%binder_history(:,this%nt1) = .false.
    
        !TODO OMP - probably not worth spawning the parallel region
        !go over all SFG units
        do m = 1, struct_group%n_elements
            !get am of m
            AM = struct_group%sfg_units(m)%get_AM(current_frame, dXdrz_db, boxdata)
            !append a m
            this%A_history(m, this%nt1) = AM(1)
            this%M_history(m, this%nt1) = AM(2)
            !update averages
            this%avA(m) = this%avA(m) + (AM(1) - this%avA(m))/this%t 
            this%avM(m) = this%avM(m) + (AM(2) - this%avM(m))/this%t 
            !subtract the best averages currently available
            this%A_history(m, this%nt1) = this%A_history(m,this%nt1) - this%avA(m)
            this%M_history(m, this%nt1) = this%M_history(m,this%nt1) - this%avM(m)

            !fill binder
            do l = 1, size(layer_selection)
                if(binder%layers(m) == layer_selection(l)) then
                    this%binder_history(m,this%nt1) = .true.
                    exit !this m was selected
                end if
            end do
        end do
    end subroutine fill_history

    subroutine skip_step(this, current_frame, struct_group, binder, dXdrz_db, layer_selection, boxdata)
        class(correlation_function_type), intent(inout) :: this
        type(current_frame_type), intent(in) :: current_frame
        type(sfg_structure_group_type), intent(in) :: struct_group
        type(group_binder_type), intent(in) :: binder
        type(dXdrz_db_type), intent(in) :: dXdrz_db
        integer(int8), dimension(:), allocatable :: layer_selection
        type(boxdata_type), intent(in) :: boxdata
        
        !now just call fill_history
        call this%fill_history(current_frame, struct_group, binder, dXdrz_db, layer_selection, boxdata)

    end subroutine skip_step

    function get_normalized(this, bd) result(res)
        class(correlation_function_type), intent(inout) :: this
        class(boxdata_type), intent(in) :: bd
        real(real64), allocatable, dimension(:) :: res
        integer :: i
        
        allocate(res(size(this%correlation_function)))
        
        do i = 1, size(this%correlation_function)
            res(i) = this%correlation_function(i) * (debye_to_ea * e_to_c) /&
                (bd%BOX_DIMENSIONS(1) * bd%BOX_DIMENSIONS(2) * 10.0_real64 * max(1,this%norm(i)))
        end do
        
    end function get_normalized
    
end module SFG_CORRELATION
