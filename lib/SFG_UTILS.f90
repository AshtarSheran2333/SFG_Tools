!------------------------------------------------------------------------------
!   SFG_UTILS
!
!   This module stores some general functions and routines that might be useful
!   throughout the whole process
!   The constants will be here too...
!
!   e.g. operator for cross product, pbc alignment...
!
!------------------------------------------------------------------------------
    
module SFG_UTILS
    use iso_fortran_env
    use BOXDATA, only: boxdata_type
    implicit none

    real(real64), parameter ::  pi = acos(-1.0_real64),& ! -
                                k_b = 1.386658e-23_real64,&
                                c = 2.99792458e10_real64,& ! cm/s
                                debye_to_ea = 0.208_real64,& ! Debye to eA
                                e_to_c = 1.6_real64 ! 1e to C

    real(real64), dimension(3), parameter ::    X_AXIS = (/1_real64, 0_real64, 0_real64/),&
                                                Y_AXIS = (/0_real64, 1_real64, 0_real64/),&
                                                Z_AXIS = (/0_real64, 0_real64, 1_real64/)
    
    complex(real64), parameter :: iunit = (0.0d0,1.0d0) 

    contains

    function pbc_minimum_image(diff, bd) result(image)
        implicit none
        real(real64), dimension(3), intent(in) :: diff
        type(boxdata_type), intent(in) :: bd
        real(real64), dimension(3) :: image

        image = diff - bd%box_dimensions * NINT( (diff) / bd%box_dimensions )

    end function pbc_minimum_image

    function pbc_wrap(vector, bd) result(back_in_the_box_baby)
        implicit none
        real(real64), dimension(3), intent(in) :: vector
        type(boxdata_type), intent(in) :: bd
        real(real64), dimension(3) :: back_in_the_box_baby

        back_in_the_box_baby = vector - bd%box_corner
        back_in_the_box_baby = modulo(back_in_the_box_baby , bd%box_dimensions) + bd%box_corner

    end function pbc_wrap

    
!TODO can be made an operator...
function cross_product(a,b)
    real(real64), dimension(3) :: cross_product
    real(real64), dimension(3), intent(IN) :: a, b

    cross_product(1) = a(2) * b(3) - a(3) * b(2)
    cross_product(2) = a(3) * b(1) - a(1) * b(3)
    cross_product(3) = a(1) * b(2) - a(2) * b(1)
end function cross_product

!dt is in fs
!omega is in cm-1
!we need to use Si units, because of unit conversions
subroutine Fourier_transform(c_t, dt, domega, omega_max, c_omega, filter)
        implicit none
        real(real64), dimension(:), intent(in) :: c_t
        real(real64), intent(in) :: dt, domega, omega_max
        complex(real64), dimension(:), allocatable, intent(inout) :: c_omega
        real(real64), intent(in), optional :: filter


        integer :: n, j, n_omega, nt
        real(real64) :: t, omega, domega_si, dt_si, weight, f, fp

        nt = size(c_t)
        n_omega = int(omega_max / domega) + 1

        if(allocated(c_omega)) deallocate(c_omega)
        allocate(c_omega(n_omega))

        c_omega = (0.0_real64, 0.0_real64)
        domega_si = 2.0_real64 * pi * c * domega ! s^-1
        dt_si = dt * 1e-15_real64 ! s
        if(present(filter)) then
            fp = filter * 1e-12 !s
        end if
        f = 1.0_real64

        !TODO can make OMP DO...
        do j = 1, n_omega

            omega = (j-1) * domega_si

            do n = 1, nt

                t = (n-1) * dt_si
                if(present(filter)) then
                    f = exp( - (t*t)/(fp*fp) )
                end if

                weight = 1.0_real64
                if(n == 1 .or. n == nt) weight = 0.5_real64

                c_omega(j) = c_omega(j) + weight * f * c_t(n) * exp(iunit * omega * t)

            end do

            c_omega(j) = -iunit * c_omega(j) * dt / omega

        end do

    end subroutine Fourier_transform

!function s_to_HHHMMSS(seconds)
!    character(14) :: s_to_HHHMMSS
!    real*8, intent(IN) :: seconds
!    integer :: hours, minutes, sec
!    
!    hours = seconds/3600
!    minutes = (seconds-hours*3600)/60
!    sec = seconds - hours*3600 - minutes*60
!    
!    write(s_to_HHHMMSS,'(I5.1,A1,I2.2,A1,I2.2)'), hours, ":", minutes, ":", sec
!    s_to_HHHMMSS = adjustl(s_to_HHHMMSS)
!end function s_to_HHHMMSS

end module
