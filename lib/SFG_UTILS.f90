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

    real(real64), parameter :: pi = acos(-1.0_real64)

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
        back_in_the_box_baby = mod(back_in_the_box_baby , bd%box_dimensions) + bd%box_corner

    end function pbc_wrap

end module