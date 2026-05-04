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
    implicit none

    real(real64), parameter :: pi = acos(-1.0_real64)

    contains
    
    function pbc_minimum_image(vector, box, corner) result(image)
        implicit none
        real(real64), dimension(3), intent(in) :: vector, box
        real(real64), dimension(3), intent(in), optional :: corner
        real(real64), dimension(3) :: image

        if(present(corner)) then 
            !shift to the box frame
            image = vector - corner
            !get minimum image
            image = image - box * ANINT( (image) / box)
            !shift the image to the correct position
            image = image + corner
        else
            !just get the minimum image
            image = vector - box * ANINT( (vector) / box)
        end if
            
    end function pbc_minimum_image
    
end module