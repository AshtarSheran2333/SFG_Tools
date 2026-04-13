module SFG_MOL_STRUCT
    use, intrinsic :: iso_fortran_env
    use FRAME_READERS, only: fr_atoms, fr_info

    !------------------the chromophore model example---------------------------
    !
    !   e.g. simple hydroxyl
    !
    !   actor -> H
    !             \
    !      base -> O
    !              |
    !              Al <- reference
    !
    !   the actors movement against base is responsible for the SFG 
    !   contribution
    !
    !   base (can be any number of atoms) serves as a reference atom for the 
    !   transformation matrix establishment average will be used in case of
    !   multiple references
    !
    !   e.g. multiple reference hydroxyl
    !
    !   actor -> H   Al
    !             \ /
    !      base -> O - () <- reference
    !               \
    !                Al
    !
    !   also note that this type just points to the actual atoms of the
    !   FRAME_READERS module
    !
    !--------------------------------------------------------------------------
    type chromophore
        integer :: actor
        integer :: base
        integer, allocatable, dimension(:) :: references
    end type
    
    !--------------------------------------------------------------------------
    !
    !   The molecule can contain any number of chromophores
    !
    !   e.g. water:
    !   H1 - actor
    !   O - base
    !   H2 - reference
    !
    !   H2 - actor
    !   O - base
    !   H1 - reference
    !
    !--------------------------------------------------------------------------
    type molecule
        type(chromophore), allocatable, dimension(:) :: chromophores
    end type
    
    contains
    
    !here probably some functions to work over some allocatable arrays returning allocated array with the molecules ???
    
    !just to test how modules work...
    subroutine pr()
        if(.not. allocated(fr_atoms)) then
            print*, "No frame in buffer"
            return
        end if
        print*, "printing from MOL_STRUCT: ", fr_atoms(1)
    end subroutine
end module