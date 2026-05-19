module BINDER_FILE
    use iso_fortran_env
    use BOXDATA
    use FRAME_READERS
    use INSTANTANEOUS_SURFACE, only: instantaneous_surface_type
    use SFG_STRUCTURE
    use SFG_UTILS
    implicit none

    integer(int64), parameter, private :: magic = Z'5245444E49424746' !BIND_SFG
    integer(int32), parameter, private :: version = 1
    
    integer(int8), dimension(8), parameter :: bLflags = (/Z'0', Z'1', Z'2', Z'3', Z'4', Z'5', Z'6', Z'7'/),&
                                            uLflags = (/Z'8', Z'9', Z'A', Z'B', Z'C', Z'D', Z'E', Z'F'/)
    
    type group_binder_type
        integer(int8), dimension(:), allocatable :: layers
    end type group_binder_type
    
    type binder_type
        type(group_binder_type), dimension(:), allocatable :: binder_groups
        character(len=32), dimension(:), allocatable :: group_names
    
    contains

    !the binder should have those methods:
    !init (from the struct)
    !write
    !read
    !verify (against struct)

    !binder frame:
    !header : magic, version, n_groups, 
    !data (groupname, n_units, binder[...])
    !data (groupname, n_units, binder[...])
    !data (groupname, n_units, binder[...])

        procedure, public :: open_file
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
    
    function open_file(this) result(res)
        class(binder_type) :: this
        integer :: res
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

    function write_frame(this) result(res)
        class(binder_type) :: this
        integer :: res
    end function

    function read_frame(this) result(res)
        class(binder_type) :: this
        integer :: res
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