module SFG_STRUCT
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
    !   The chromophore can have none of the references, then reference for calculations
    !   must be somehow preselected, e.g. Z axis
    !
    !       base -> C - H <- actor
    !
    !   the chromophore is the smallest "unit of spectrum" that we can get
    !   the chromophore should carry information about what set of parameters will be used for the A - M calculation
    !
    !--------------------------------------------------------------------------
    type chromophore
        integer :: actor
        integer :: base
        integer :: parameters_id
        integer, allocatable, dimension(:) :: references
    end type chromophore
    
    !--------------------------------------------------------------------------
    !
    !   The SFG unit can contain any number of chromophores
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
    type sfg_unit
        type(chromophore), allocatable, dimension(:) :: chromophores
    end type sfg_unit
    
    !----------------------------the input-------------------------------------
    !
    !   the input should look like this:
    !
    !   chromophore:
    !   PAR_ID  BASE    ACTOR   REFERENCE(S) - up to 10
    !   ID      ID      ID      ID  ... ID
    !
    !   a special group for waters ($WATERS):
    !   PAR_ID  O   H1  H2
    !   ID      ID  ID  ID
    !
    !   a special group for hydroxyls ($HYDROXYLS):
    !   chromophore
    !   chromophore
    !
    !   a special group for other ($OTHER)
    !   - where user can specify groups as a set of chromophores
    !   - the group can have a name, the groups with the same name will be stored in a separate list
    !   $GROUP          -
    !   chromophore     |
    !   chromophore     |---- default group
    !   chromophore     |
    !   $GROUP          -
    !   chromophore     ----- simple chromophore (default group)
    !   $GROUP CH3      -
    !   chromophore     |
    !   chromophore     |---- CH3 group
    !   chromophore     |
    !   $GROUP CH3      -
    !   chromophore     ----- simple chromophore (default group)
    !   chromophore     ----- simple chromophore (default group)
    !   chromophore     ----- simple chromophore (default group)
    !   $GROUP CH3      -
    !   chromophore     |
    !   chromophore     |---- CH3 group
    !   chromophore     |
    !   $GROUP          -
    !   $GROUP MONSTER  -
    !   chromophore     |
    !   chromophore     |
    !   chromophore     |
    !   chromophore     |---- MONSTER group
    !   chromophore     |
    !   chromophore     |
    !   chromophore     |
    !   $GROUP          -
    !
    !   this should yield a structure:
    !
    !   | DEFAULT | CH3 | MONSTER |
    !   |    5    |  2  |    1    |
    !   ---------------------------
    !
    !                   the file can look like this:
    !   $WATERS
    !   ID ID ID ID
    !   ID ID ID ID
    !   ID ID ID ID
    !
    !   $HYDROXYLS
    !   ID ID ID ID ... ID
    !   ID ID ID ID ... ID
    !   ID ID ID ID ... ID
    !
    !   $OTHER
    !   ID ID ID ID ... ID
    !   $GROUP
    !   ID ID ID ID ... ID
    !   ID ID ID ID ... ID
    !   ID ID ID ID ... ID
    !   ID ID ID ID ... ID
    !   $GROUP
    !   ID ID ID ID ... ID
    !   
    !--------------------------------------------------------------------------
    
    type structgroup
        integer(kind = int32) :: n_elements
        character(len = 32) :: name
        type(sfg_unit), dimension(:), allocatable :: sfg_units
    end type
    
    integer, parameter, private ::  STATE_GROUP = 0,&
                                    STATE_NONE = 1,&
                                    STATE_WATERS = 2,&
                                    STATE_HYDROXYLS = 3,&
                                    STATE_OTHER = 4
    
    integer, private :: file
    integer, private :: line_number
    
    type(structgroup), protected, dimension(:), allocatable :: sfg_structure
    integer, protected :: sfg_structure_max_index

    private :: clear_struct
    private :: read_water
    private :: clear_chgroup_container
    private :: append_chgroup

    
    contains

    !! this function returns maximal index, this can be used to check agains frame_reader n_group to be safe when using those together
    !! returns maximal index stored in the water, hydroxyl, and other array
    !function sfg_group_struct_get_max_index() result(res)
    !    integer :: res
    !    integer :: i,j,k
    !    
    !    !this is pretty ugly and slow, but whatever
    !    res = -1

    !    !water
    !    do i = 1, n_water
    !        do j = 1, size(water(i)%chromophores)
    !            if(res < water(i)%chromophores(j)%base) res = water(i)%chromophores(j)%base
    !            if(res < water(i)%chromophores(j)%actor) res = water(i)%chromophores(j)%actor
    !            do k = 1, size(water(i)%chromophores(j)%references)
    !                if(res < water(i)%chromophores(j)%references(k)) res = water(i)%chromophores(j)%references(k)
    !            end do
    !        end do
    !    end do
    !    
    !    !hydroxyl
    !    do i = 1, n_hydroxyl
    !        do j = 1, size(hydroxyl(i)%chromophores)
    !            if(res < hydroxyl(i)%chromophores(j)%base) res = hydroxyl(i)%chromophores(j)%base
    !            if(res < hydroxyl(i)%chromophores(j)%actor) res = hydroxyl(i)%chromophores(j)%actor
    !            do k = 1, size(hydroxyl(i)%chromophores(j)%references)
    !                if(res < hydroxyl(i)%chromophores(j)%references(k)) res = hydroxyl(i)%chromophores(j)%references(k)
    !            end do
    !        end do
    !    end do

    !    !other
    !    do i = 1, n_other
    !        do j = 1, size(other(i)%chromophores)
    !            if(res < other(i)%chromophores(j)%base) res = other(i)%chromophores(j)%base
    !            if(res < other(i)%chromophores(j)%actor) res = other(i)%chromophores(j)%actor
    !            do k = 1, size(other(i)%chromophores(j)%references)
    !                if(res < other(i)%chromophores(j)%references(k)) res = other(i)%chromophores(j)%references(k)
    !            end do
    !        end do
    !    end do

    !end function sfg_group_struct_get_max_index
    !
    !subroutine clear_struct()
    !    if(allocated(water)) deallocate(water)
    !    n_water = 0
    !    if(allocated(hydroxyl)) deallocate(hydroxyl)
    !    n_hydroxyl = 0
    !    if(allocated(other)) deallocate(other)
    !    n_other = 0
    !end subroutine clear_struct
    !
    !subroutine clear_group_container(dest, dest_length)
    !    type(chgroup), allocatable, dimension(:), intent(inout) :: dest
    !    integer, intent(inout) :: dest_length

    !    dest_length = 0

    !    if(.not. allocated(dest)) return
    !    
    !    deallocate(dest)
    !end subroutine clear_group_container
    !
    
    !append a group with name groupname
    !if that group does not exist, create that group
    !todo make private
    subroutine append_group(groupname, unit)
        implicit none
        character(*), intent(in) :: groupname
        type(sfg_unit), intent(in) :: unit

        type(structgroup), allocatable, dimension(:) :: temp_group
        type(sfg_unit), allocatable, dimension(:) :: temp_units

        integer :: i, newsize
        
        do i=1, size(sfg_structure)
            if(sfg_structure(i)%name == trim(adjustl(groupname))) exit !got the group
        end do

        if(i > size(sfg_structure)) then
        !need to alloc new structgroup
            newsize = size(sfg_structure) + 1
            allocate(temp_group(newsize))
            temp_group(1:size(sfg_structure)) = sfg_structure
            call move_alloc(from=temp_group, to=sfg_structure)
        end if

        !just push to struct i 
        if(.not. allocated(sfg_structure(i)%sfg_units)) then
            allocate(sfg_structure(i)%sfg_units(1))
            sfg_structure(i)%name = trim(adjustl(groupname))
            sfg_structure(i)%n_elements = 1
            sfg_structure(i)%sfg_units = unit
        else
            newsize = size(sfg_structure(i)%sfg_units) + 1
            allocate(temp_units(newsize))
            temp_units(1:size(sfg_structure(i)%sfg_units)) = sfg_structure(i)%sfg_units
            temp_units(newsize:newsize) = unit
            call move_alloc(from = temp_units, to=sfg_structure(i)%sfg_units)
            sfg_structure(i)%n_elements = sfg_structure(i)%n_elements + 1
        end if
    end subroutine append_group
    
    ! result = 0 - OK
    ! result = -1 - incomplete group
    ! result = -2 - IO error
    function read_water(group) result(res)
        type(sfg_unit), intent(inout) :: group
        integer :: res
        integer, dimension(4) :: ids
        character(128) :: line
        logical :: is_open
        
        res = -2
        inquire(file, opened = is_open)
        if(.not. is_open) return
        
        if(allocated(group%chromophores)) deallocate(group%chromophores)
        allocate(group%chromophores(2))
        allocate(group%chromophores(1)%references(1))
        allocate(group%chromophores(2)%references(1))

        read(file, "(A)", iostat = ierr) line
        res = -2; if(ierr .ne. 0) return
        line_number = line_number + 1
        
        read(line, *, iostat = ierr) ids
        if(ierr .ne. 0) then !incomplete group
            res = -1
            backspace(file) !this can be either another section, or wrong line
            line_number = line_number - 1
            return
        end if

        group%chromophores(1)%parameters_id = ids(1)
        group%chromophores(2)%parameters_id = ids(1)
        
        group%chromophores(1)%base = ids(2)
        group%chromophores(2)%base = ids(2)
        
        group%chromophores(1)%actor = ids(3)
        group%chromophores(2)%actor = ids(4)

        group%chromophores(1)%references(1) = ids(4)
        group%chromophores(2)%references(1) = ids(3)

        if(maxval(ids(2:)) > sfg_structure_max_index) sfg_structure_max_index = maxval(ids(2:))
        res = 0
    end function read_water

    ! result = 0 - OK
    ! result = -1 - incomplete group
    ! result = -2 - IO error
    function read_simple_chromophore_chgroup(group) result(res)
        type(sfg_unit), intent(inout) :: group
        integer :: res
        integer, parameter :: max_hydroxyl_references = 12
        integer, dimension(max_hydroxyl_references) :: ids
        integer :: i, ierr
        character(128) :: line
        logical :: is_open
        
        res = -2
        inquire(file, opened = is_open)
        if(.not. is_open) return
        
        if(allocated(group%chromophores)) deallocate(group%chromophores)
        allocate(group%chromophores(1))

        read(file, "(A)", iostat = ierr) line
        res = -2; if(ierr .ne. 0) return
        line_number = line_number + 1

        read(line, *, iostat = ierr) (ids(i), i = 1, max_hydroxyl_references)
        if(i <= 3) then !incomplete chromophore - must have PAR_ID, ACTOR, BASE
            backspace(file)
            line_number = line_number - 1
            res = -1
            return
        end if
        

        group%chromophores(1)%parameters_id = ids(1)

        group%chromophores(1)%base = ids(2)
        
        group%chromophores(1)%actor = ids(3)

        !any references?
        i = i - 1 !i was incremented by extra one to either finish the loop or when read failed
        if( (i-3) .gt. 0 ) then 
            allocate(group%chromophores(1)%references(i-3))
            group%chromophores(1)%references = ids(4:i)
        end if
        
        if(maxval(ids(2:i)) > sfg_structure_max_index) sfg_structure_max_index = maxval(ids(2:i))
        res = 0
    end function read_simple_chromophore_chgroup

    subroutine append_chgroup_chromophores(group, chr)
        type(sfg_unit), intent(inout) :: group
        type(chromophore), intent(in) :: chr
        type(chromophore), dimension(:), allocatable :: temp_chromophores
        integer :: group_chromophores_size
        
        if(.not. allocated(group%chromophores)) then
            allocate(group%chromophores(1))
            group%chromophores(1) = chr
            return
        end if
        
        group_chromophores_size = size(group%chromophores)

        allocate(temp_chromophores(group_chromophores_size + 1))
        temp_chromophores(1:group_chromophores_size) = group%chromophores
        temp_chromophores(group_chromophores_size + 1) = chr
        
        deallocate(group%chromophores)
        call move_alloc(from=temp_chromophores, to=group%chromophores)
    end subroutine append_chgroup_chromophores
    
    function read_other(group) result(res)
        type(sfg_unit), intent(inout) :: group
        integer :: res
        logical :: is_open
        type(sfg_unit) :: temp_group
        
        res = -2
        inquire(file, opened = is_open)
        if(.not. is_open) return
        
        if(allocated(group%chromophores)) then
            res = read_simple_chromophore_chgroup(temp_group)
            if(res .ne. 0) return
            call append_chgroup_chromophores(group, temp_group%chromophores(1))
        else
            !read simple chromophore
            res = read_simple_chromophore_chgroup(group)
            if(res .ne. 0) return
        end if
    end function read_other

    
    ! result STATE_NONE - STATE_OTHER - OK
    ! result = -2 - IO error
    function get_reader_state(groupname) result(res)
        integer :: res
        character(len=32), intent(inout) :: groupname
        character(128) :: line
        character(64) :: W1
        integer :: pos
        
        read(file, "(A)", iostat = ierr) line 
        res = -2; if(ierr .ne. 0) return 
        line_number = line_number + 1
        
        line = trim(adjustl(line))
        
        pos = index(line, ' ')
        if(pos .gt. 0) then
            W1 = line(1:pos)
            groupname = line(pos+1:len(line))
        else
            W1 = line
            groupname = ""
        end if

        select case(W1)
            case ("$WATERS")
                res = STATE_WATERS    
            case ("$HYDROXYLS")
                res = STATE_HYDROXYLS
            case ("$OTHER")
                res = STATE_OTHER
            case ("$GROUP")
                res = STATE_GROUP
                if( (trim(adjustl(groupname)) == "$WATERS")&
                    .or. (trim(adjustl(groupname)) == "$HYDROXYLS")&
                    .or. (trim(adjustl(groupname)) == "$WATERS") ) then
                    
                    write(error_unit, "(A,' ',A,' (line: ',I0,')')") "WARNING: SFG_STRUCT invalid groupname $GROUP", trim(groupname), line_number
                    write(error_unit, "(A)") "Assigned groupname: $UNASSIGNED"
                    groupname = "$UNASSIGNED"
                end if
            case default
                if( (len_trim(adjustl(W1)) .gt. 0) ) then !not an empty line
                    if( index(trim(adjustl(W1)), '#') .ne. 1 ) then !does not start with #
                        write(error_unit, "(A,' ',A,' (line: ',I0,')')") "WARNING: SFG_STRUCT discarded:", trim(line), line_number
                    end if
                end if
                res = STATE_NONE
        end select
    end function get_reader_state

    ! result = 0 - OK
    ! result = -1 - wrong format
    ! result = -2 - IO error
    function read_struct(filename) result(res)
        !this is disgusting...
        character(*) :: filename
        integer :: res
        logical :: is_open
        integer :: read_state = STATE_NONE
        type(sfg_unit) :: group
        integer :: ret
        character(128) :: line
        logical :: reading_group = .false.
        character(32) :: groupname, current_groupname

        inquire(file, opened = is_open)
        if(is_open) close(file)
        line_number = 0
        sfg_structure_max_index = -1
        
        open(newunit = file, file = filename, status = 'old', iostat = ierr)
        res = -2; if(ierr .ne. 0) return

        if(allocated(sfg_structure)) deallocate(sfg_structure)
        
        !get reader to some state...
        do while(read_state <= STATE_NONE)
            read_state = get_reader_state(groupname)
            res = -2; if(read_state < STATE_NONE) return !did not find any label
        end do
        
        do while(read_state > STATE_NONE)
            select case(read_state)
                case (STATE_WATERS)
                    ret = read_water(group)    
                    if(ret .eq. -2) exit !IO error - get out of the loop
                    if(ret .eq. 0) then !append group and continue reading
                        call append_group("$WATERS", group)
                        cycle
                    end if
                    read_state = get_reader_state(groupname)
                    if(read_state < 0) exit !IO error - get out of the loop
                    if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                        read_state = STATE_WATERS
                        cycle
                    end if
                    deallocate(group%chromophores) !switch state
                case (STATE_HYDROXYLS)
                    ret = read_simple_chromophore_chgroup(group)    
                    if(ret .eq. -2) exit !IO error - get out of the loop
                    if(ret .eq. 0) then !append group and continue reading
                        call append_group("$HYDROXYLS", group)
                        cycle
                    end if
                    read_state = get_reader_state(groupname)
                    if(read_state < 0) exit !IO error - get out of the loop
                    if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                        read_state = STATE_HYDROXYLS
                        cycle
                    end if
                    deallocate(group%chromophores) !switch state
                case (STATE_OTHER)
                    if(.not. reading_group) then
                        ret = read_other(group)    
                        if(ret .eq. -2) exit !IO error - get out of the loop
                        if(ret .eq. 0) then !append group and continue reading
                            !simple chromophore without group
                            call append_group("$OTHER", group)
                            if(allocated(group%chromophores)) deallocate(group%chromophores)
                            cycle
                        end if
                        read_state = get_reader_state(groupname)
                        if(read_state < 0) exit !IO error - get out of the loop
                        if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                            if(read_state == STATE_GROUP) then
                                current_groupname = groupname
                                reading_group = .not. reading_group
                            end if
                            read_state = STATE_OTHER
                            if(allocated(group%chromophores)) deallocate(group%chromophores)
                            cycle
                        end if
                        if(allocated(group%chromophores)) deallocate(group%chromophores) !another state
                    else
                        ret = read_other(group)    
                        if(ret .eq. -2) exit !IO error - get out of the loop
                        if(ret .eq. 0) cycle !chromophores appended 
                        read_state = get_reader_state(groupname)
                        if(read_state < 0) exit !IO error - get out of the loop
                        if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                            if(read_state == STATE_GROUP) then
                                !some logic with current group name...
                                !if group with groupname exists, append, else, 
                                if(len_trim(current_groupname) == 0) current_groupname = "$OTHER"
                                call append_group(current_groupname, group)
                                reading_group = .not. reading_group
                                if(allocated(group%chromophores)) deallocate(group%chromophores)
                            end if
                            read_state = STATE_OTHER
                            cycle
                        end if
                        if(reading_group) then
                            res = -1 !wrong format
                            return
                        end if
                        if(allocated(group%chromophores)) deallocate(group%chromophores) !another state
                    end if
                case default
                    error_stop("ERROR: read_struct state machine error")
            end select
        end do
        
        res = -2; if((n_water + n_hydroxyl + n_other) .eq. 0) return !didnt get any data from the file
        res = 0
    end function read_struct
    
    !here probably some functions to work over some allocatable arrays returning allocated array with the groups ???
    
end module
