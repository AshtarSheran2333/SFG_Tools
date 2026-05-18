module BINDER_FILE
    use iso_fortran_env
    implicit none

    integer(int64), parameter :: magic = Z'5245444E49424746' !BIND_SFG
    integer(int32), parameter :: version = 1

    type binder_type
        integer(int8), dimension(:,:), allocatable :: binder
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
    end type binder_type

    contains
    
    function open_file(this) result(res)
        class(binder_type) :: this
        integer :: res
    end function

    function init(this) result(res)
        class(binder_type) :: this
        integer :: res
    end function

    function write_frame(this) result(res)
        class(binder_type) :: this
        integer :: res
    end function

    function read_frame(this) result(res)
        class(binder_type) :: this
        integer :: res
    end function

    function verify_struct(this) result(res)
        class(binder_type) :: this
        integer :: res
    end function
    
end module