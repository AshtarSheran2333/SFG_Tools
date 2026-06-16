module SFG_STRUCTURE
    use, intrinsic :: iso_fortran_env
    use FRAME_READERS, only: current_frame_type
    use DXDRZ_DB, only: dXdrz_db_type
    use BOXDATA, only: boxdata_type
    use SFG_UTILS, only: pbc_minimum_image, cross_product
    use UTILS_ERROR

#include "utils_error_macros.h"

    !------------------the site model example---------------------------
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
    !   The site can have none of the references, then reference for calculations
    !   must be somehow preselected, e.g. Z axis
    !
    !       base -> C - H <- actor
    !
    !   the site is the smallest "unit of spectrum" that we can get
    !   the site should carry information about what set of parameters will be used for the A - M calculation
    !
    !--------------------------------------------------------------------------
    type sfg_site_type
        integer :: actor
        integer :: base
        integer :: parameters_id
        integer, allocatable, dimension(:) :: references
    end type sfg_site_type
    
    !--------------------------------------------------------------------------
    !
    !   The SFG unit can contain any number of sites
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
    type sfg_unit_type
        type(sfg_site_type), allocatable, dimension(:) :: sites
        integer(kind = int32) :: n_unique_bases
        integer, dimension(:), allocatable :: unique_bases !holds unique bases -> simpler density calculations...
    
    contains
    
        procedure, public :: fill_unique_bases
        procedure, public :: append_site
        procedure, public :: get_AM
    end type sfg_unit_type
    
    !----------------------------the input-------------------------------------
    !
    !   the input should look like this:
    !
    !   site:
    !   PAR_ID  BASE    ACTOR   REFERENCE(S) - up to 10
    !   ID      ID      ID      ID  ... ID
    !
    !   a special group for waters ($WATERS):
    !   PAR_ID  O   H1  H2
    !   ID      ID  ID  ID
    !
    !   a special group for hydroxyls ($HYDROXYLS):
    !   site
    !   site
    !
    !   a special group for other ($OTHER)
    !   - where user can specify groups as a set of sites
    !   - the group can have a name, the groups with the same name will be stored in a separate list
    !   $GROUP          -
    !   site     |
    !   site     |---- default group
    !   site     |
    !   $GROUP          -
    !   site     ----- simple site (default group)
    !   $GROUP CH3      -
    !   site     |
    !   site     |---- CH3 group
    !   site     |
    !   $GROUP CH3      -
    !   site     ----- simple site (default group)
    !   site     ----- simple site (default group)
    !   site     ----- simple site (default group)
    !   $GROUP CH3      -
    !   site     |
    !   site     |---- CH3 group
    !   site     |
    !   $GROUP          -
    !   $GROUP MONSTER  -
    !   site     |
    !   site     |
    !   site     |
    !   site     |---- MONSTER group
    !   site     |
    !   site     |
    !   site     |
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
    
    type sfg_structure_group_type
        integer(kind = int32) :: n_elements
        character(len = 32) :: name
        type(sfg_unit_type), dimension(:), allocatable :: sfg_units
    end type sfg_structure_group_type
    
    integer, parameter, private ::  STATE_GROUP = 0,&
                                    STATE_NONE = 1,&
                                    STATE_WATERS = 2,&
                                    STATE_HYDROXYLS = 3,&
                                    STATE_OTHER = 4
    
    type sfg_structure_type
        !TODO this should be protected, not to be overwritten from outside - no time to write getters...
        type(sfg_structure_group_type), dimension(:), allocatable :: groups
        integer :: max_atom_index

        integer, private :: file
        integer, private :: line_number

    contains
        procedure, public :: read_structure

        procedure, private :: clear_structure
        procedure, private :: append_group
        procedure, private :: read_water
        procedure, private :: read_simple_site
        procedure, private :: read_other

        procedure, private :: get_reader_state
        
    end type sfg_structure_type

    
    contains

    subroutine clear_structure(this)
        implicit none
        class(sfg_structure_type), intent(inout) :: this

        if(allocated(this%groups)) deallocate(this%groups)
    end subroutine clear_structure
    
    subroutine fill_unique_bases(this)
        class(sfg_unit_type), intent(inout) :: this
        integer, dimension(:), allocatable :: temp_unique_bases
        integer :: i,j
        logical :: exists

        if(allocated(this%unique_bases)) then
            deallocate(this%unique_bases)
            this%n_unique_bases = 0
        end if
        
        this%n_unique_bases = 0

        !go through all the bases of SFG unit
        do i = 1, size(this%sites)
            if(.not. allocated(this%unique_bases)) then !first iteration
                allocate(this%unique_bases(1))
                this%unique_bases(1) = this%sites(i)%base
                this%n_unique_bases = 1
                cycle
            end if

            !is this base unique?
            exists = any(this%unique_bases == this%sites(i)%base) 
            if(exists) cycle
            
            !found unique base - append unique_bases
            allocate(temp_unique_bases(this%n_unique_bases + 1))
            
            temp_unique_bases(1:this%n_unique_bases) = this%unique_bases
            temp_unique_bases(this%n_unique_bases + 1) = this%sites(i)%base
            
            call move_alloc(from = temp_unique_bases, to = this%unique_bases)
            this%n_unique_bases = this%n_unique_bases + 1
        end do
        
    end subroutine fill_unique_bases
    
    
    !append a group with name groupname
    !if that group does not exist, create that group
    !makes sure that the unique bases are filled
    subroutine append_group(this, groupname, sfg_unit)
        implicit none
        class(sfg_structure_type), intent(inout) :: this
        character(*), intent(in) :: groupname
        type(sfg_unit_type), intent(inout) :: sfg_unit

        type(sfg_structure_group_type), allocatable, dimension(:) :: temp_group
        type(sfg_unit_type), allocatable, dimension(:) :: temp_units

        integer :: i, newsize
        
        do i=1, size(this%groups)
            if(this%groups(i)%name == trim(adjustl(groupname))) exit !got the group
        end do

        if(i > size(this%groups)) then
        !need to alloc new group
            newsize = size(this%groups) + 1
            allocate(temp_group(newsize))
            temp_group(1:size(this%groups)) = this%groups
            call move_alloc(from=temp_group, to=this%groups)
        end if

        call sfg_unit%fill_unique_bases()

        !just push to struct i 
        if(.not. allocated(this%groups(i)%sfg_units)) then
            allocate(this%groups(i)%sfg_units(1))
            this%groups(i)%name = trim(adjustl(groupname))
            this%groups(i)%n_elements = 1
            this%groups(i)%sfg_units = sfg_unit
        else
            newsize = size(this%groups(i)%sfg_units) + 1
            allocate(temp_units(newsize))
            temp_units(1:size(this%groups(i)%sfg_units)) = this%groups(i)%sfg_units
            temp_units(newsize:newsize) = sfg_unit
            call move_alloc(from = temp_units, to=this%groups(i)%sfg_units)
            this%groups(i)%n_elements = this%groups(i)%n_elements + 1
        end if
    end subroutine append_group
    
    ! result = 0 - OK
    ! result = -1 - incomplete group
    ! result = -2 - IO error
    function read_water(this, sfg_unit) result(res)
        implicit none
        class(sfg_structure_type), intent(inout) :: this
        type(sfg_unit_type), intent(inout) :: sfg_unit
        integer :: res
        integer, dimension(4) :: ids
        character(128) :: line
        logical :: is_open
        integer :: ierr
        
        res = -2
        inquire(this%file, opened = is_open)
        if(.not. is_open) return
        
        if(allocated(sfg_unit%sites)) deallocate(sfg_unit%sites)
        allocate(sfg_unit%sites(2))
        allocate(sfg_unit%sites(1)%references(1))
        allocate(sfg_unit%sites(2)%references(1))

        read(this%file, "(A)", iostat = ierr) line
        res = -2; if(ierr .ne. 0) return
        this%line_number = this%line_number + 1
        
        read(line, *, iostat = ierr) ids
        if(ierr .ne. 0) then !incomplete sfg_unit
            res = -1
            backspace(this%file) !this can be either another section, or wrong line
            this%line_number = this%line_number - 1
            return
        end if

        sfg_unit%sites(1)%parameters_id = ids(1)
        sfg_unit%sites(2)%parameters_id = ids(1)
        
        sfg_unit%sites(1)%base = ids(2)
        sfg_unit%sites(2)%base = ids(2)
        
        sfg_unit%sites(1)%actor = ids(3)
        sfg_unit%sites(2)%actor = ids(4)

        sfg_unit%sites(1)%references(1) = ids(4)
        sfg_unit%sites(2)%references(1) = ids(3)

        if(maxval(ids(2:)) > this%max_atom_index) this%max_atom_index = maxval(ids(2:))
        res = 0
    end function read_water

    ! result = 0 - OK
    ! result = -1 - incomplete group
    ! result = -2 - IO error
    function read_simple_site(this, sfg_unit) result(res)
        implicit none
        class(sfg_structure_type), intent(inout) :: this
        type(sfg_unit_type), intent(inout) :: sfg_unit
        integer :: res
        integer, parameter :: max_site_references = 12
        integer, dimension(max_site_references) :: ids
        integer :: i, ierr
        character(128) :: line
        logical :: is_open
        
        res = -2
        inquire(this%file, opened = is_open)
        if(.not. is_open) return
        
        if(allocated(sfg_unit%sites)) deallocate(sfg_unit%sites)
        allocate(sfg_unit%sites(1))

        read(this%file, "(A)", iostat = ierr) line
        res = -2; if(ierr .ne. 0) return
        this%line_number = this%line_number + 1

        ids = 0
        read(line, *, iostat = ierr) (ids(i), i = 1, max_site_references)
        if(count(ids .ne. 0) < 3) then !incomplete site - must have PAR_ID, ACTOR, BASE
            backspace(this%file)
            this%line_number = this%line_number - 1
            res = -1
            return
        end if

        sfg_unit%sites(1)%parameters_id = ids(1)

        sfg_unit%sites(1)%base = ids(2)
        
        sfg_unit%sites(1)%actor = ids(3)

        !any references?
        i = i - 1 !i was incremented by extra one to either finish the loop or when read failed
        if( (i-3) .gt. 0 ) then 
            allocate(sfg_unit%sites(1)%references(i-3))
            sfg_unit%sites(1)%references = ids(4:i)
        end if
        
        if(maxval(ids(2:i)) > this%max_atom_index) this%max_atom_index = maxval(ids(2:i))
        res = 0
    end function read_simple_site

    subroutine append_site(this, sfg_site)
        implicit none
        class(sfg_unit_type), intent(inout) :: this
        type(sfg_site_type), intent(in) :: sfg_site
        type(sfg_site_type), dimension(:), allocatable :: temp_sites
        integer :: n_sites
        
        if(.not. allocated(this%sites)) then
            allocate(this%sites(1))
            this%sites(1) = sfg_site
            return
        end if
        
        n_sites = size(this%sites)

        allocate(temp_sites(n_sites + 1))
        temp_sites(1:n_sites) = this%sites
        temp_sites(n_sites + 1) = sfg_site
        
        deallocate(this%sites)
        call move_alloc(from=temp_sites, to=this%sites)
    end subroutine append_site
    
    !gets A and M of SFG unit
    function get_AM(this, frame, param_db, boxdata) result(res)
        implicit none
        class(sfg_unit_type), intent(in) :: this
        type(current_frame_type), intent(in) :: frame
        type(dXdrz_db_type), intent(in) :: param_db
        type(boxdata_type), intent(in) :: boxdata
        real(real64), dimension(2) :: res

        integer :: n, ref, i, j
        real(real64), dimension(3) :: u, v, diff!base-actor, base-reference
        real(real64), dimension(3,3) :: D
        real(real64) :: r, scale, A, M, vz

        res = 0
        A = 0
        M = 0
        
        !loop over the sites
        do n = 1, size(this%sites)    
            if( (this%sites(n)%parameters_id > param_db%count) .or. &
                (this%sites(n)%parameters_id > param_db%count) ) &
                error_stop("invalid parameters_id") !TODO this should be checked somewhere else...

            !CONSTRUCT D MATRIX
            !solve actor - base vector
            u = frame%positions(:,this%sites(n)%actor) - frame%positions(:,this%sites(n)%base)
            u = pbc_minimum_image(u, boxdata)
            r = norm2(u)
            ! z component
            D(:,3) = u / r

            !solve base-reference vector
            i = size(this%sites(n)%references)
            if(i == 0) then
                !no base -> use -Z
                v = (/0.0, 0.0, -1.0/)
            else
                !average the bases
                v = 0
                do ref = 1, i
                    diff = frame%positions(:,this%sites(n)%references(ref)) - frame%positions(:,this%sites(n)%base)
                    diff = pbc_minimum_image(diff, boxdata)
                    v = v + diff
                end do
                v = v / i
            end if

            ! x component
            scale = dot_product(v,D(:,3))
            D(:,1) = scale*D(:,3)-v
            r = norm2(D(:,1))
            if(r .ne. 0) then 
                !everything OK
                D(:,1) = D(:,1) / r
            else
                !troubles with division by 0
                if((D(1,3) == 0) .and. (D(2,3) == 0)) then
                    ! be careful, in this case we have no clue about the direction
                    ! leads to division by zero
                    ! X axis if lab Z == Z else -X axis
                    D(:,1) = (/D(3,3),0d0,0d0/)
                else
                    ! if actor base reference angle is 180deg
                    ! chose X perpendicular to Z
                    D(:,1) = (/D(2,3),-D(1,3),0d0/)
                    r = norm2(D(:,1))
                    D(:,1) = D(:,1) / r
                end if
            end if
            
            ! y component
            D(:,2) = cross_product(D(:,3),D(:,1))
            !TODO end construct D matrix subroutine

            !get vz
            diff = frame%velocities(:,this%sites(n)%actor) - frame%velocities(:,this%sites(n)%base)
            vz = dot_product(diff,D(:,3))

            !TODO calculate the A M
            !M_R
            do i=1,3     ! do on x,y,z
                associate( dMdrz => param_db%dMdrz_record(this%sites(n)%parameters_id)%elements(i) )
                M = M + ( D(boxdata%R,i) * dMdrz * vz )
                end associate
            end do
        
            !A_PQ
            do i=1,3     ! do on x,y,z
                do j=1,3     ! do on x,y,z
                    associate( dAdrz => param_db%dAdrz_record(this%sites(n)%parameters_id)%elements(i,j) )
                    A = A +( D(boxdata%P,i) * dAdrz * D(boxdata%Q,j) * vz )
                    end associate
                end do
            end do

        end do

        res = (/A, M/)

    end function get_AM

    function read_other(this, sfg_unit) result(res)
        implicit none
        class(sfg_structure_type), intent(inout) :: this
        type(sfg_unit_type), intent(inout) :: sfg_unit
        integer :: res
        logical :: is_open
        type(sfg_unit_type) :: temp_sfg_unit
        
        res = -2
        inquire(this%file, opened = is_open)
        if(.not. is_open) return
        
        if(allocated(sfg_unit%sites)) then
            res = this%read_simple_site(temp_sfg_unit)
            if(res .ne. 0) return
            call sfg_unit%append_site(temp_sfg_unit%sites(1))
        else
            !read simple site
            res = this%read_simple_site(sfg_unit)
            if(res .ne. 0) return
        end if
    end function read_other

    
    ! result STATE_NONE - STATE_OTHER - OK
    ! sets the groupname if specified $TOKEN GROUPNAME
    ! result = -2 - IO error
    function get_reader_state(this, groupname) result(res)
        implicit none
        class(sfg_structure_type), intent(inout) :: this
        integer :: res
        character(len=32), intent(inout) :: groupname
        character(128) :: line
        character(64) :: W1
        integer :: pos
        integer :: ierr
        
        read(this%file, "(A)", iostat = ierr) line 
        res = -2; if(ierr .ne. 0) return 
        this%line_number = this%line_number + 1
        
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
                    
                    write(error_unit, "(A,' ',A,' (line: ',I0,')')") "WARNING: SFG_STRUCT invalid groupname $GROUP",&
                        trim(groupname), this%line_number
                    write(error_unit, "(A)") "Assigned groupname: $UNASSIGNED"
                    groupname = "$UNASSIGNED"
                end if
            case default
                if( (len_trim(adjustl(W1)) .gt. 0) ) then !not an empty line
                    if( index(trim(adjustl(W1)), '#') .ne. 1 ) then !does not start with #
                        write(error_unit, "(A,' ',A,' (line: ',I0,')')") "WARNING: SFG_STRUCT discarded:", trim(line), this%line_number
                    end if
                end if
                res = STATE_NONE
        end select
    end function get_reader_state

    ! result = 0 - OK
    ! result = -1 - wrong format
    ! result = -2 - IO error
    function read_structure(this, filename) result(res)
        implicit none
        !this is disgusting... the whole file parsing is crammed here... whatever
        class(sfg_structure_type), intent(inout) :: this
        character(*) :: filename
        integer :: res
        logical :: is_open
        integer :: read_state = STATE_NONE
        type(sfg_unit_type) :: sfg_unit
        integer :: ret
        character(128) :: line
        logical :: reading_group = .false.
        character(32) :: groupname, current_groupname
        integer :: ierr

        inquire(this%file, opened = is_open)
        if(is_open) close(this%file)
        this%line_number = 0
        this%max_atom_index = -1
        
        open(newunit = this%file, file = filename, status = 'old', form = 'formatted', access = 'sequential' , iostat = ierr)
        res = -2; if(ierr .ne. 0) return

        call this%clear_structure()
        
        !get reader to some state...
        do while(read_state <= STATE_NONE)
            read_state = this%get_reader_state(groupname)
            res = -2; if(read_state < STATE_NONE) return !did not find any label
        end do
        
        do while(read_state > STATE_NONE)
            select case(read_state)
                case (STATE_WATERS)
                    ret = this%read_water(sfg_unit)    
                    if(ret .eq. -2) exit !IO error - get out of the loop
                    if(ret .eq. 0) then !append group and continue reading
                        call this%append_group("$WATERS", sfg_unit)
                        cycle
                    end if
                    read_state = this%get_reader_state(groupname)
                    if(read_state < 0) exit !IO error - get out of the loop
                    if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                        read_state = STATE_WATERS
                        cycle
                    end if
                    deallocate(sfg_unit%sites) !switch state
                case (STATE_HYDROXYLS)
                    ret = this%read_simple_site(sfg_unit)    
                    if(ret .eq. -2) exit !IO error - get out of the loop
                    if(ret .eq. 0) then !append group and continue reading
                        call this%append_group("$HYDROXYLS", sfg_unit)
                        cycle
                    end if
                    read_state = this%get_reader_state(groupname)
                    if(read_state < 0) exit !IO error - get out of the loop
                    if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                        read_state = STATE_HYDROXYLS
                        cycle
                    end if
                    deallocate(sfg_unit%sites) !switch state
                case (STATE_OTHER)
                    if(.not. reading_group) then
                        ret = this%read_other(sfg_unit)    
                        if(ret .eq. -2) exit !IO error - get out of the loop
                        if(ret .eq. 0) then !append group and continue reading
                            !simple site without group
                            call this%append_group("$OTHER", sfg_unit)
                            if(allocated(sfg_unit%sites)) deallocate(sfg_unit%sites)
                            cycle
                        end if
                        read_state = this%get_reader_state(groupname)
                        if(read_state < 0) exit !IO error - get out of the loop
                        if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                            if(read_state == STATE_GROUP) then
                                current_groupname = groupname
                                reading_group = .not. reading_group
                            end if
                            read_state = STATE_OTHER
                            if(allocated(sfg_unit%sites)) deallocate(sfg_unit%sites)
                            cycle
                        end if
                        if(allocated(sfg_unit%sites)) deallocate(sfg_unit%sites) !another state
                    else
                        ret = this%read_other(sfg_unit)    
                        if(ret .eq. -2) exit !IO error - get out of the loop
                        if(ret .eq. 0) cycle !sites appended 
                        read_state = this%get_reader_state(groupname)
                        if(read_state < 0) exit !IO error - get out of the loop
                        if(read_state <= STATE_NONE) then !just a garbage line, continue reading
                            if(read_state == STATE_GROUP) then
                                !some logic with current group name...
                                !if group with groupname exists, append, else, 
                                if(len_trim(current_groupname) == 0) current_groupname = "$OTHER"
                                call this%append_group(current_groupname, sfg_unit)
                                reading_group = .not. reading_group
                                if(allocated(sfg_unit%sites)) deallocate(sfg_unit%sites)
                            end if
                            read_state = STATE_OTHER
                            cycle
                        end if
                        if(reading_group) then
                            res = -1 !wrong format
                            return
                        end if
                        if(allocated(sfg_unit%sites)) deallocate(sfg_unit%sites) !another state
                    end if
                case default
                    error_stop("ERROR: read_struct state machine error")
            end select
        end do
        
        res = -2; if(.not. allocated(this%groups)) return !didnt get any data from the file
        res = 0
    end function read_structure
    
    !here probably some functions to work over some allocatable arrays returning allocated array with the groups ???
    
end module
