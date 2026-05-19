module BINDER_FILE
    use iso_fortran_env
    use BOXDATA
    use FRAME_READERS
    use INSTANTANEOUS_SURFACE, only: instantaneous_surface_type
    use SFG_STRUCTURE
    use SFG_UTILS
    implicit none

    integer(int64), parameter, private :: binder_magic = Z'4746535F444E4942' !BIND_SFG
    integer(int32), parameter, private :: binder_version = 1
    
    integer(int8), dimension(8), parameter :: bLflags = (/Z'0', Z'1', Z'2', Z'3', Z'4', Z'5', Z'6', Z'7'/),&
                                            uLflags = (/Z'8', Z'9', Z'A', Z'B', Z'C', Z'D', Z'E', Z'F'/)
    
    type group_binder_type
        integer(int8), dimension(:), allocatable :: layers
    end type group_binder_type
    
    type binder_type
        type(group_binder_type), dimension(:), allocatable :: binder_groups
        character(len=32), dimension(:), allocatable :: group_names

        integer, private :: file_unit
    
    contains

    !the binder should have those methods:
    !init (from the struct)
    !write
    !read
    !verify (against struct)

        procedure, public :: open_file
        procedure, public :: close_file
        procedure, public :: init
        procedure, public :: write_frame
        procedure, public :: read_frame
        procedure, public :: verify_struct
        procedure, public :: fill_binder
    end type binder_type

    contains
    
    !return -1 struct does not match the internal binder structure
    !return 0 - OK
    function fill_binder(this, fr, struct, instasurf, bd, z_center) result(res)
        implicit none
        class(binder_type), intent(inout) :: this
        class(frame_reader), intent(in) :: fr
        class(sfg_structure_type), intent(in) :: struct
        class(instantaneous_surface_type), intent(in) :: instasurf
        real(real64), intent(in) :: z_center
        class(boxdata_type), intent(in) :: bd
        integer :: res

        integer :: group, unit, base, layer, i
        real(real64), dimension(3) :: position
        real(real64), dimension(2) :: distances
        real(real64) :: distance

        res = 0
        
        !if not allocated, try to allocate
        if(.not. allocated(this%binder_groups)) res = this%init(struct)
        if(res .ne. 0) return
        
        !verify the dimensions
        res = this%verify_struct(struct)
        if(res .ne. 0) return
        
        !TODO can be made parallel, but is it worth it?
        do group = 1, size(struct%groups)
            do unit = 1, struct%groups(group)%n_elements
                position = 0
                do base = 1, struct%groups(group)%sfg_units(unit)%n_unique_bases
                    associate( base_id => struct%groups(group)%sfg_units(unit)%unique_bases(base) )
                    position = position + fr%frame%positions(:, base_id)
                    end associate
                end do
                position = position / struct%groups(group)%sfg_units(unit)%n_unique_bases
                position = pbc_wrap(position, bd)
                
                distances = instasurf%get_distances(position, bd)
                
                !get distance, based on Z look at 1 or 2
                !assign layer, save it
                if(position(3) <= z_center) then
                    distance = distances(1)
                    layer = size(bd%LAYERS_LIMITS) + 1
                    do i = 1, size(bd%LAYERS_LIMITS)
                        if(distance < bd%LAYERS_LIMITS(i)) then
                            layer = i
                            exit
                        end if
                    end do
                    this%binder_groups(group)%layers(unit) = bLFlags(layer)
                else
                    distance = distances(2)
                    layer = size(bd%LAYERS_LIMITS) + 1
                    do i = 1, size(bd%LAYERS_LIMITS)
                        if(distance < bd%LAYERS_LIMITS(i)) then
                            layer = i
                            exit
                        end if
                    end do
                    this%binder_groups(group)%layers(unit) = uLFlags(layer)
                end if
            
            end do
        end do
    
    end function
    
    function open_file(this, read_only, must_exist, name) result(res)
        class(binder_type), intent(inout) :: this
        character(*), intent(in), optional :: name
        logical, intent(in), optional :: read_only, must_exist
        character(128) :: file_name
        integer :: res
        logical :: is_open
        character(20) :: act, stat
        
        inquire(this%file_unit, opened = is_open)
        if(is_open) close(this%file_unit)
        
        file_name = "binder.bin"
        stat = 'UNKNOWN'
        act = 'READWRITE'
        
        if(present(name)) then
            file_name = trim(adjustl(name))
        end if
        
        if(present(must_exist)) then
            if(must_exist) stat = 'OLD'
        end if

        if(present(read_only)) then
            if(read_only) act = 'READ'
        end if
            
        open(newunit = this%file_unit, &
                file = trim(adjustl(file_name)), &
                status = stat, &
                action = act, &
                form = "unformatted", &
                access = "stream", &
                convert = "little_endian", &
                iostat = res)
    end function

    function close_file(this) result(res)
        class(binder_type), intent(inout) :: this
        integer :: res
        logical :: is_open
        
        inquire(this%file_unit, opened = is_open)
        if(is_open) close(this%file_unit)
        res = 0
    end function

    function init(this, struct) result(res)
        class(binder_type) :: this
        class(sfg_structure_type), intent(in) :: struct
        integer :: res
        integer :: group
        
        if(allocated(this%binder_groups)) deallocate(this%binder_groups)
        if(allocated(this%group_names)) deallocate(this%group_names)

        allocate(this%binder_groups(size(struct%groups)))
        allocate(this%group_names(size(struct%groups)))
        do group = 1, size(struct%groups)
            allocate(this%binder_groups(group)%layers(struct%groups(group)%n_elements))
            this%binder_groups(group)%layers = Z'FF' !illegal
        end do

        this%group_names = struct%groups(:)%name
        
        res = 0 !allways return 0 not expectiong fail in the allocate - which is naive I guess
    end function

    !binder frame:
    !header : magic(int64), version(int32), n_groups(int32), 
    !(n_units(int32), groupname(32))[n_groups]
    !data(each nibble = 1 record - same order)[n_groups] - in case of odd records, append 0 nibble
    function write_frame(this) result(res)
        class(binder_type) :: this
        logical :: is_open
        integer :: res, group, unit
        integer(int8) :: chunk
        
        inquire(this%file_unit, opened = is_open)
        res = -1; if(.not. is_open) return

        res = 0

        !write header
        if(res == 0) write(this%file_unit, iostat = res) binder_magic
        if(res == 0) write(this%file_unit, iostat = res) binder_version
        if(res == 0) write(this%file_unit, iostat = res) int(size(this%binder_groups),int32)
        do group = 1, size(this%binder_groups)
            if(res == 0) write(this%file_unit, iostat = res) int(size(this%binder_groups(group)%layers),int32)
            if(res == 0) write(this%file_unit, iostat = res) this%group_names(group)
        end do
        
        !write data
        do group = 1, size(this%binder_groups)
            do unit = 1, size(this%binder_groups(group)%layers)
                if(mod(unit,2) == 1) then
                    chunk = 0
                    chunk = SHIFTL(this%binder_groups(group)%layers(unit),4)
                else
                    chunk = OR(chunk, AND(this%binder_groups(group)%layers(unit), Z'F'))
                    if(res == 0) write(this%file_unit, iostat = res) chunk
                end if
            end do
            !in case of odd number of units, write the last one - just lowest nibble
            if(mod(unit,2) == 0) then
                if(res == 0) write(this%file_unit, iostat = res) chunk
            end if
        end do

    end function
    
    !TODO expecting that the binder will be initiated from struct...
    !TODO make the binder initialize just from the file...
    !result -1 - file not open
    !result -2 - unexpected format
    !result < 0 - IO-error
    !result 0 - OK
    function read_frame(this) result(res)
        class(binder_type) :: this
        integer :: res
        integer(int64) :: magic
        integer(int32) :: i32, groups, group_size, i, n_nibbles, j
        integer(int8) :: i8
        character(len=32) :: groupname
        logical :: is_open, is_odd
        
        inquire(this%file_unit, opened = is_open)
        res = -1; if(.not. is_open) return
        
        read(this%file_unit, iostat = res) magic !read magic
        if(res .ne. 0) return
        if(magic .ne. binder_magic) then
            res = -2
            return
        end if

        read(this%file_unit, iostat = res) i32 !read version
        if(res .ne. 0) return
        if(i32 .ne. binder_version) then
            res = -2
            return
        end if

        read(this%file_unit, iostat = res) groups !read n_groups
        if(res .ne. 0) return
        if( (groups .ne. size(this%binder_groups)) .or. &
            (groups .ne. size(this%group_names)) ) then
            res = -2
            return
        end if

        do i = 1, groups
            read(this%file_unit, iostat = res) i32 !read group_size
            if(res .ne. 0) return
            if(i32 .ne. size(this%binder_groups(i)%layers)) then
                res = -2
                return
            end if
            
            read(this%file_unit, iostat = res) groupname !read group_name
            if(res .ne. 0) return
            if(groupname .ne. this%group_names(i)) then
                res = -2
                return
            end if
        end do

        !header OK

        !TODO THIS IS NOT OK... I AM TIRED
        do i = 1, groups
            n_nibbles = ceiling(real(size(this%binder_groups(i)%layers)) / 2.0) 
            is_odd = .false.
            if(mod(size(this%binder_groups(i)%layers), 2) == 1) is_odd = .true.
            
            do j = 1, n_nibbles
                read(this%file_unit, iostat = res) i8 !read nibbles
                !read lower nibble
                this%binder_groups(i)%layers(j*2-1) = OR(0, AND(i8, Z'F'))   
                if(j == n_nibbles .and. is_odd) exit
                !read upper nibble
                this%binder_groups(i)%layers(j*2) = OR(0, SHIFTR(i8, 4))
            end do
        end do
        print*, "ok"

            
        
        
        !if(.not. allocated(this%binder_groups))
    end function

    function verify_struct(this, struct) result(res)
        class(binder_type) :: this
        class(sfg_structure_type), intent(in) :: struct
        integer :: res, group

        res = -1
        !verify struct
        if(size(this%binder_groups) .ne. size(struct%groups)) return
        if(any(this%group_names .ne. struct%groups(:)%name)) return
        do group = 1, size(struct%groups)
            if(size(this%binder_groups(group)%layers) .ne. struct%groups(group)%n_elements) return    
        end do

        res = 0
    end function
    
end module