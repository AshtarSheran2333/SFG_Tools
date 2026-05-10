module DENSITY
    use iso_fortran_env
    use UTILS_ERROR
    implicit none
#include "utils_error_macros.h"
        
    !again, all members public, makes my life easier - do not touch them
    type density_profile_type
        integer(int64) :: n_frames
        real(real64) :: x_min, x_max, delta, bin_volume
        real(real64), dimension(:), allocatable :: bins
        character(64) :: name
        
    contains
        procedure, public :: init
        procedure, public :: add_point
        procedure, public :: next_frame
        procedure, public :: write_file
       
    end type

    contains
    
    subroutine init(this, x_min, x_max, n_bins, box_dimensions, name)
        implicit none
        class(density_profile_type), intent(inout) :: this
        real(real64), intent(in) :: x_min, x_max
        real(real64), dimension(3), intent(in) :: box_dimensions
        integer(int64), intent(in) :: n_bins
        integer(int64) :: nb
        character(*), intent(in), optional :: name
        !save the parameters

        if(x_min == x_max) error_stop("x_min == x_max") !problem
        
        if(x_min < x_max) then
            this%x_min = x_min
            this%x_max = x_max
        else
            this%x_min = x_max
            this%x_max = x_min
        end if 
        
        nb = n_bins
        if(nb <= 0) nb = 1 !in case of wrong input

        this%delta = (this%x_max - this%x_min) / real(nb, real64)
        this%bin_volume = box_dimensions(1) * box_dimensions(2) * this%delta

        this%n_frames = 0
        
        if(present(name)) then
            this%name = name
        else
            this%name = ""
        end if
        
        !allocate the array based on the parameters
        if(allocated(this%bins)) deallocate(this%bins)
        allocate(this%bins(nb))

        this%bins = 0.0_real64
    end subroutine

    subroutine add_point(this, x,y)
        implicit none
        class(density_profile_type), intent(inout) :: this
        real(real64), intent(in) :: x,y
        integer :: pos

        !assign the bins position
        pos = int((x - this%x_min) / this%delta) + 1
        
        if(pos < 1 .or. pos > size(this%bins)) return !does not fit, exit
        
        !add to the bins
        !!$omp atomic update - this can be done to make all members private, but it might affect performance
        this%bins(pos) = this%bins(pos) + y
    end subroutine

    subroutine next_frame(this)
        implicit none
        class(density_profile_type), intent(inout) :: this

        this%n_frames = this%n_frames + 1
    end subroutine

    subroutine write_file(this, name)
        implicit none
        character(64), optional :: name !overrides this%name if present
        class(density_profile_type), intent(inout) :: this
        integer :: file, i
        !write file, do not forget to divide each bin by number of steps
        !if name "" do not write the file...
        open(newunit = file, file = trim(adjustl(this%name))//"_density.dat")
        
        do i = 1, size(this%bins)
            write(file, *) real(i-1)*this%delta + this%x_min, this%bins(i)/(this%n_frames * this%bin_volume)
        end do
        
        close(file)
        !TODO not finished
    end subroutine

end module DENSITY