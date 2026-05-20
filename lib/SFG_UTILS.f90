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
                                k_b = 1.386658d-23,&
                                c = 2.99792458d10,& ! cm/s
                                debye_to_ea = 0.208,& ! Debye to eA
                                e_to_c = 1.6 ! 1e to C

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

function cross_product(a,b)
	real(real64), dimension(3) :: cross_product
	real(real64), dimension(3), intent(IN) :: a, b

	cross_product(1) = a(2) * b(3) - a(3) * b(2)
	cross_product(2) = a(3) * b(1) - a(1) * b(3)
	cross_product(3) = a(1) * b(2) - a(2) * b(1)
end function cross_product

!subroutine fourier integral...

!function s_to_HHHMMSS(seconds)
!	character(14) :: s_to_HHHMMSS
!	real*8, intent(IN) :: seconds
!	integer :: hours, minutes, sec
!	
!	hours = seconds/3600
!	minutes = (seconds-hours*3600)/60
!	sec = seconds - hours*3600 - minutes*60
!	
!	write(s_to_HHHMMSS,'(I5.1,A1,I2.2,A1,I2.2)'), hours, ":", minutes, ":", sec
!	s_to_HHHMMSS = adjustl(s_to_HHHMMSS)
!end function s_to_HHHMMSS

end module