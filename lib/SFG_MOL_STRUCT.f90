module SFG_MOL_STRUCT
    use, intrinsic :: iso_fortran_env
    use FRAME_READERS, only: fr_atoms, fr_info
    use UTILS_ERROR

#include "utils_error_macros.h"

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
    end type chromophore
    
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
    end type molecule
    
    !----------------------------the input-------------------------------------
    !
    !   the input should look like this:
    !
    !   chromophore:
    !   BASE    ACTOR   REFERENCE(S)
    !   ID      ID      ID  ... ID
    !
    !   a special group for waters:
    !   O   H1  H2
    !   ID  ID  ID
    !
    !   a special group for hydroxyls:
    !   chromophore
    !   chromophore
    !
    !   a special group for others
    !   - where user can specify molecules as a set of chromophores
    !   MOL             |
    !   chromophore     |
    !   chromophore     |---- molecule
    !   chromophore     |
    !   MOL             |
    !   chromophore     ----- simple chromophore
    !
    !                   the file can look like this:
    !   $WATERS
    !   ID ID ID
    !   ID ID ID
    !   ID ID ID
    !
    !   $HYDROXYLS
    !   ID ID ID ... ID
    !   ID ID ID ... ID
    !   ID ID ID ... ID
    !
    !   $OTHER
    !   ID ID ID ... ID
    !   MOL
    !   ID ID ID ... ID
    !   ID ID ID ... ID
    !   ID ID ID ... ID
    !   ID ID ID ... ID
    !   MOL
    !   ID ID ID ... ID
    !--------------------------------------------------------------------------
    
    integer, parameter, private ::  STATE_MOL = 0,&
                                    STATE_NONE = 1,&
                                    STATE_WATERS = 2,&
                                    STATE_HYDROXYLS = 3,&
                                    STATE_OTHER = 4
    
    integer, private :: file
    
    integer, private, parameter :: alloc_chunk = 100
    
    integer, protected :: n_water = 0
    type(molecule), allocatable, dimension(:), protected :: water
    integer, protected :: n_hydroxyl = 0
    type(molecule), allocatable, dimension(:), protected :: hydroxyl
    integer, protected :: n_other = 0
    type(molecule), allocatable, dimension(:), protected :: other

    private :: clear_struct
    private :: read_water
    private :: clear_molecule_container
    private :: append_molecule

    
    contains
    
    subroutine clear_struct()
        if(allocated(water)) deallocate(water)
        n_water = 0
        if(allocated(hydroxyl)) deallocate(hydroxyl)
        n_hydroxyl = 0
        if(allocated(other)) deallocate(other)
        n_other = 0
    end subroutine clear_struct
    
    subroutine clear_molecule_container(dest, dest_length)
        type(molecule), allocatable, dimension(:), intent(inout) :: dest
        integer, intent(inout) :: dest_length

        dest_length = 0

        if(.not. allocated(dest)) return
        
        deallocate(dest)
    end subroutine clear_molecule_container
    
    subroutine append_molecule(dest, dest_length, mol)
        implicit none
        type(molecule), allocatable, dimension(:), intent(inout) :: dest
        integer, intent(inout) :: dest_length
        type(molecule), intent(in) :: mol
        type(molecule), allocatable, dimension(:) :: temp
        integer :: new_len
        
        if(.not. allocated(dest)) then
            allocate(dest(alloc_chunk))
            dest_length = 0
        end if

        if(dest_length >= size(dest)) then !reallocate
            new_len = size(dest) + alloc_chunk
            allocate(temp(new_len))

            temp(1:dest_length) = dest
            
            call move_alloc(from=temp, to=dest)
        end if

        dest_length = dest_length + 1
        dest(dest_length) = mol
    end subroutine append_molecule
    
    ! result = 0 - OK
    ! result = -1 - incomplete molecule
    ! result = -2 - IO error
    function read_water(mol) result(res)
        type(molecule), intent(inout) :: mol
        integer :: res
        integer :: ID1, ID2, ID3
        character(128) :: line
        logical :: is_open
        
        res = -2
        inquire(file, opened = is_open)
        if(.not. is_open) return
        
        if(allocated(mol%chromophores)) deallocate(mol%chromophores)
        allocate(mol%chromophores(2))
        allocate(mol%chromophores(1)%references(1))
        allocate(mol%chromophores(2)%references(1))

        read(file, "(A)", iostat = ierr) line
        res = -2; if(ierr .ne. 0) return
        
        read(line, *, iostat = ierr) ID1, ID2, ID3
        if(ierr .ne. 0) then !incomplete molecule
            res = -1
            backspace(file) !this can be either another section, or wrong line
            return
        end if

        mol%chromophores(1)%base = ID1
        mol%chromophores(2)%base = ID1
        
        mol%chromophores(1)%actor = ID2
        mol%chromophores(2)%actor = ID3

        mol%chromophores(1)%references(1) = ID3
        mol%chromophores(2)%references(1) = ID2
        res = 0
    end function read_water

    ! result = 0 - OK
    ! result = -1 - incomplete molecule
    ! result = -2 - IO error
    function read_simple_chromophore_mol(mol) result(res)
        type(molecule), intent(inout) :: mol
        integer :: res
        integer, parameter :: max_hydroxyl_references = 12
        integer, dimension(max_hydroxyl_references) :: ids
        integer :: i, ierr
        character(128) :: line
        logical :: is_open
        
        res = -2
        inquire(file, opened = is_open)
        if(.not. is_open) return
        
        if(allocated(mol%chromophores)) deallocate(mol%chromophores)
        allocate(mol%chromophores(1))

        read(file, "(A)", iostat = ierr) line
        res = -2; if(ierr .ne. 0) return

        read(line, *, iostat = ierr) (ids(i), i = 1, max_hydroxyl_references)
        if(i <= 3) then !incomplete hydroxyl
            backspace(file)
            res = -1
            return
        end if
        
        i = i - 1
        allocate(mol%chromophores(1)%references(i-2))
        mol%chromophores(1)%base = ids(1)
        
        mol%chromophores(1)%actor = ids(2)

        mol%chromophores(1)%references = ids(3:i)

        res = 0
    end function read_simple_chromophore_mol

    subroutine append_mol_chromophores(mol, chr)
        type(molecule), intent(inout) :: mol
        type(chromophore), intent(in) :: chr
        type(chromophore), dimension(:), allocatable :: temp_chromophores
        integer :: mol_chromophores_size
        
        if(.not. allocated(mol%chromophores)) then
            allocate(mol%chromophores(1))
            mol%chromophores(1) = chr
            return
        end if
        
        mol_chromophores_size = size(mol%chromophores)

        allocate(temp_chromophores(mol_chromophores_size + 1))
        temp_chromophores(1:mol_chromophores_size) = mol%chromophores
        temp_chromophores(mol_chromophores_size + 1) = chr
        
        deallocate(mol%chromophores)
        call move_alloc(from=temp_chromophores, to=mol%chromophores)
    end subroutine append_mol_chromophores
    
    function read_other(mol) result(res)
        type(molecule), intent(inout) :: mol
        integer :: res
        integer, parameter :: max_hydroxyl_references = 12
        integer, dimension(max_hydroxyl_references) :: ids
        integer :: i, ierr
        character(128) :: line
        logical :: is_open
        type(molecule) :: temp_mol
        
        res = -2
        inquire(file, opened = is_open)
        if(.not. is_open) return
        
        if(allocated(mol%chromophores)) then
            res = read_simple_chromophore_mol(temp_mol)
            if(res .ne. 0) return
            call append_mol_chromophores(mol, temp_mol%chromophores(1))
        else
            !read simple chromophore
            res = read_simple_chromophore_mol(mol)
            if(res .ne. 0) return
        end if
    end function read_other

    
    ! result STATE_NONE - STATE_OTHER - OK
    ! result = -2 - IO error
    function get_reader_state() result(res)
        integer :: res
        character(128) :: line

        read(file, "(A)", iostat = ierr) line 
        res = -2; if(ierr .ne. 0) return 

        line = trim(line)
        
        select case(line)
            case ("$WATERS")
                res = STATE_WATERS    
            case ("$HYDROXYLS")
                res = STATE_HYDROXYLS
            case ("$OTHER")
                res = STATE_OTHER
            case ("MOL")
                res = STATE_MOL
            case default
                res = STATE_NONE
        end select
    end function get_reader_state

    ! result = 0 - OK
    ! result = -1 - wrong format
    ! result = -2 - IO error
    function read_struct(filename) result(res)
        character(*) :: filename
        integer :: res
        logical :: is_open
        integer :: read_state = STATE_NONE
        type(molecule) :: mol
        integer :: ret
        character(128) :: line
        logical :: reading_mol = .false.

        inquire(file, opened = is_open)
        if(is_open) close(file)
        
        open(newunit = file, file = filename, status = 'old', iostat = ierr)
        res = -2; if(ierr .ne. 0) return

        call clear_struct()
        
        !get reader to some state...
        do while(read_state <= STATE_NONE)
            read_state = get_reader_state()
            res = -2; if(read_state < STATE_NONE) return !did not find any label
        end do
        
        do while(read_state > STATE_NONE)
            select case(read_state)
                case (STATE_WATERS)
                    ret = read_water(mol)    
                    if(ret .eq. -2) exit !IO error - get out of the loop
                    if(ret .eq. 0) then !append molecule and continue reading
                        call append_molecule(water, n_water, mol)
                        cycle
                    end if
                    read_state = get_reader_state()
                    if(read_state < 0) exit !IO error - get out of the loop
                    if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                        read_state = STATE_WATERS
                        cycle
                    end if
                    deallocate(mol%chromophores) !switch state
                case (STATE_HYDROXYLS)
                    ret = read_simple_chromophore_mol(mol)    
                    if(ret .eq. -2) exit !IO error - get out of the loop
                    if(ret .eq. 0) then !append molecule and continue reading
                        call append_molecule(hydroxyl, n_hydroxyl, mol)
                        cycle
                    end if
                    read_state = get_reader_state()
                    if(read_state < 0) exit !IO error - get out of the loop
                    if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                        read_state = STATE_HYDROXYLS
                        cycle
                    end if
                    deallocate(mol%chromophores) !switch state
                case (STATE_OTHER)
                    if(.not. reading_mol) then
                        ret = read_other(mol)    
                        if(ret .eq. -2) exit !IO error - get out of the loop
                        if(ret .eq. 0) then !append molecule and continue reading
                            call append_molecule(other, n_other, mol)
                            cycle
                        end if
                        read_state = get_reader_state()
                        if(read_state < 0) exit !IO error - get out of the loop
                        if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                            if(read_state == STATE_MOL) then
                                reading_mol = .not. reading_mol
                            end if
                            read_state = STATE_OTHER
                            if(allocated(mol%chromophores)) deallocate(mol%chromophores)
                            cycle
                        end if
                        if(allocated(mol%chromophores)) deallocate(mol%chromophores) !another state
                    else
                        ret = read_other(mol)    
                        if(ret .eq. -2) exit !IO error - get out of the loop
                        if(ret .eq. 0) cycle !chromophores appended 
                        read_state = get_reader_state()
                        if(read_state < 0) exit !IO error - get out of the loop
                        if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                            if(read_state == STATE_MOL) then
                                call append_molecule(other, n_other, mol)
                                reading_mol = .not. reading_mol
                                if(allocated(mol%chromophores)) deallocate(mol%chromophores)
                            end if
                            read_state = STATE_OTHER
                            cycle
                        end if
                        if(reading_mol) then
                            res = -1 !wrong format
                            return
                        end if
                        if(allocated(mol%chromophores)) deallocate(mol%chromophores) !another state
                    end if
                case default
                    error_stop("ERROR: read_struct state machine error")
            end select
        end do
        
        res = -2; if((n_water + n_hydroxyl + n_other) .eq. 0) return !didnt get any data from the file
        res = 0
    end function read_struct
    
    !here probably some functions to work over some allocatable arrays returning allocated array with the molecules ???
    
    !just to test how modules work...
    subroutine pr()
        type(molecule), allocatable :: pushed
        integer :: i

        i = read_struct("struct")
        print*, "read struct res: ", i
    end subroutine
end module
