module INSTANTANEOUS_SURFACE
    use iso_fortran_env
    use FRAME_READERS, only: current_frame_type
    use SFG_UTILS, only: pbc_minimum_image, pbc_wrap, pi
    use BOXDATA, only: boxdata_type
    use UTILS_ERROR
    implicit none
    
#include "utils_error_macros.h"
    
    type :: instantaneous_surface_type
        integer(int64), dimension(3) :: n_points !number of iterations through space
        real(real64), dimension(3) :: volume_element
        real(real64), dimension(3) :: start !corner of the mesh
        real(real64), dimension(:,:), allocatable :: up_mesh
        real(real64), dimension(:,:), allocatable :: bot_mesh
        integer(int64), dimension(:,:), allocatable :: up_index !index of iteration to get the point
        integer(int64), dimension(:,:), allocatable :: bot_index !index of iteration to get the point

        integer, private :: instasurf_bin_unit, instasurf_xyz_unit

    contains
    
        procedure, public :: init
        procedure, public :: init_flat
        procedure, public :: calculate
        procedure, public :: get_distances !result (/bot distance, up distance/)
        procedure, public :: get_center
        procedure, public :: write_grid_interface
        procedure, public :: open_bin_file
        procedure, public :: open_xyz_file
        procedure, public :: write_bin_frame
        procedure, public :: write_xyz_frame
        procedure, public :: read_next
        procedure, public :: close_bin_file
        procedure, public :: close_xyz_file
        
    end type instantaneous_surface_type
    
    integer(int32), parameter :: instantaneous_surface_magic = INT(Z'53474653', int32) !will be "SFGS" when hexdumped
    integer(int32), parameter :: instantaneous_surface_bin_version = 0
    
contains

    subroutine init(this, bd)
        implicit none
        class(instantaneous_surface_type), intent(inout) :: this
        type(boxdata_type), intent(in) :: bd
        
        !wipe
        if(allocated(this%up_mesh)) deallocate(this%up_mesh)
        if(allocated(this%bot_mesh)) deallocate(this%bot_mesh)
        if(allocated(this%up_index)) deallocate(this%up_index)
        if(allocated(this%bot_index)) deallocate(this%bot_index)
        
        !get mesh
        this%n_points = NINT(bd%box_dimensions/bd%interface_volume_element)

        this%volume_element = bd%box_dimensions/this%n_points

        this%start = bd%box_corner - this%volume_element
        
        !allocate containers
        
        associate( i => this%n_points(1), &
                   j => this%n_points(2))

        allocate(this%up_mesh(i,j))
        allocate(this%bot_mesh(i,j))
        allocate(this%up_index(i,j))
        allocate(this%bot_index(i,j))

        end associate

        this%up_mesh = 0.0
        this%bot_mesh = 0.0
        this%up_index = this%n_points(3)
        this%bot_index = 1
        
    end subroutine init

    subroutine init_flat(this, bd, bot, up)
        implicit none
        class(instantaneous_surface_type), intent(inout) :: this
        type(boxdata_type), intent(in) :: bd
        real(real64), intent(in) :: bot, up

        call this%init(bd)

        this%up_mesh = up
        this%bot_mesh = bot
    end subroutine init_flat

    function calculate(this, frame, bd) result(res)
        implicit none
        !TODO - ugly implementatation - calling the same block twice
        class(instantaneous_surface_type), intent(inout) :: this
        type(current_frame_type), intent(in) :: frame
        type(boxdata_type), intent(in) :: bd
        integer(int32) :: res
        real(real64) :: density_threshold, tollerance

        integer(int64) :: i,j,k,m
        logical :: found_up_interface
        logical :: found_bot_interface
        real(real64), dimension(3) :: pos, prev_pos, diff
        real(real64) :: rho, rhodiff, prev_rhodiff, r

        density_threshold = bd%liquid_bulk_number_density/2.0
        tollerance = 0.004 !bd%liquid_bulk_number_density * 0.1 - could be like this

        res = 0
        
        !$OMP PARALLEL DO &
        !$OMP DEFAULT(NONE) &
        !$OMP SHARED(this, frame, bd, density_threshold, tollerance) &
        !$OMP PRIVATE(i, j, k, m, found_up_interface, found_bot_interface, pos, prev_pos, diff, rho, rhodiff, prev_rhodiff, r) &
        !$OMP REDUCTION(+:res)
        do i=1,this%n_points(1)
            res = 0
            do j=1,this%n_points(2)
            
                found_up_interface = .false.
                found_bot_interface = .false.
                
                pos = this%start
                rhodiff = density_threshold
    
                !probing from bottom
                do k=max(1_int64,this%bot_index(i,j)-bd%interface_pushback), this%n_points(3)
                    prev_pos = pos
                    pos = this%volume_element * (/i,j,k/) + this%start
            
                    !evaluate rho
                    rho = 0
            
                    do m=1, size(bd%liquid_heavy_atoms)
                        associate( atom_pos => frame%positions(:,bd%liquid_heavy_atoms(m)) )
                        diff = atom_pos - pos
                        end associate
                        !pbc correction              
                        diff = pbc_minimum_image(diff, bd)
                        r = norm2(diff)
                        !cutoff after 3 sigma, the value would be too small, save some calculation time
                        if( r <= 3*bd%coarse_graining_length ) then
                            rho = rho + exp(-r**2/(2*bd%coarse_graining_length**2))/((2*pi*bd%coarse_graining_length**2)**1.5)
                        end if  
                    end do
            
                    !shaping the rho
                    prev_rhodiff = rhodiff
                    rhodiff = abs(density_threshold - rho)
                    
                    !if the positive derivative is found -> we have found the interface
                    if(prev_rhodiff < tollerance) then
                        if(rhodiff > prev_rhodiff) then
                            found_bot_interface = .true.
                            this%bot_mesh(i,j) = prev_pos(3)
                            this%bot_index(i,j) = k - 1
                            exit
                        end if
                    end if
                end do
                
                pos = this%start
                rhodiff = density_threshold
    
                !poking rod from top
                do k=min(this%n_points(3), this%up_index(i,j)+bd%interface_pushback),1,-1
                    prev_pos = pos
                    pos = this%volume_element * (/i,j,k/) + this%start
            
                    !evaluate rho
                    rho = 0
            
                    do m=1, size(bd%liquid_heavy_atoms)
                        associate( atom_pos => frame%positions(:,bd%liquid_heavy_atoms(m)) )
                        diff = atom_pos - pos
                        end associate
                        !pbc correction              
                        diff = pbc_minimum_image(diff, bd)
                        r = norm2(diff)
                        !cutoff after 3 sigma, the value would be too small, save some calculation time
                        if( r <= 3*bd%coarse_graining_length ) then
                            rho = rho + exp(-r**2/(2*bd%coarse_graining_length**2))/((2*pi*bd%coarse_graining_length**2)**1.5)
                        end if  
                    end do

                    !shaping the rho
                    prev_rhodiff = rhodiff
                    rhodiff = abs(density_threshold - rho)
            
                    !if the positive derivative is found -> we have found the interface
                    if(prev_rhodiff < tollerance) then
                        if(rhodiff > prev_rhodiff) then
                        found_up_interface = .true.
                        this%up_mesh(i,j) = prev_pos(3)
                        this%up_index(i,j) = k - 1
                        exit
                        end if
                    end if

                end do
                
                !checks if program found both interfaces
        
                if( .not. found_up_interface) res = res - 1
                if( .not. found_bot_interface) res = res - 1
            end do
        end do
        !$OMP END PARALLEL DO

        if(res < 0) write(error_unit,"(A)") "ERROR: instasurf did not find all the points of interface grid"
        !TODO - here we potentially can look for the neighbor points and try to interpolate the missing points, but this should not happen
    end function calculate

    function get_distances(this, point, bd) result(distances)
        implicit none
        class(instantaneous_surface_type), intent(in) :: this
        real(real64), dimension(3), intent(in) :: point
        type(boxdata_type), intent(in) :: bd
        real(real64), dimension(3) :: pos
        real(real64) :: xgrid, ygrid, xt, yu
        real(real64), dimension(2) :: distances
        integer(int32) :: x, x1, y, y1
        
        pos = pbc_wrap(point, bd)
        
        !find where the point belongs on the grid
        !X
        x = int((pos(1) - (this%start(1) + this%volume_element(1))) / this%volume_element(1)) + 1
        x1 = x + 1
        !Y
        y = int((pos(2) - (this%start(2) + this%volume_element(2))) / this%volume_element(2)) + 1
        y1 = y + 1

        !deal with the edge cases
        if(x > this%n_points(1)) x = 1
        if(x1 > this%n_points(1)) x1 = 1
        if(y > this%n_points(2)) y = 1
        if(y1 > this%n_points(2)) y1 = 1
        
        xgrid = this%start(1) + x*this%volume_element(1)
        ygrid = this%start(2) + y*this%volume_element(2)
        
        xt = (pos(1)-xgrid)/(this%volume_element(1))
        yu = (pos(2)-ygrid)/(this%volume_element(2))
        
        !interface z
        distances(2) = (1-xt)*(1-yu)*this%up_mesh(x,y) &
                    + xt*(1-yu)*this%up_mesh(x1,y) &
                    + (1-xt)*yu*this%up_mesh(x,y1) &
                    + xt*yu*this%up_mesh(x1,y1)

        distances(1) = (1-xt)*(1-yu)*this%bot_mesh(x,y) &
                    + xt*(1-yu)*this%bot_mesh(x1,y) &
                    + (1-xt)*yu*this%bot_mesh(x,y1) &
                    + xt*yu*this%bot_mesh(x1,y1)
        !distance from the interface
        distances(1) = pos(3) - distances(1)
        distances(2) = distances(2) - pos(3)
    end function get_distances

    function get_center(this) result(res)
        class(instantaneous_surface_type), intent(in) :: this
        real(real64) :: res

        if( (size(this%bot_mesh) .eq. 0) .or. (size(this%up_mesh) .eq. 0)) then
            error_stop("empty instantaneous surface - cannot get a center")
        end if

        res = sum(this%bot_mesh) + sum(this%up_mesh)
        res = res / (size(this%bot_mesh) + size(this%up_mesh))
    end function get_center

    function write_grid_interface(this, filename) result(res)
        class(instantaneous_surface_type), intent(inout) :: this
        character(*), intent(IN), optional :: filename
        integer :: ierr, file_unit, res
        logical :: is_open
        
        if(present(filename)) then
            open(newunit = file_unit, file=trim(adjustl(filename)), iostat = ierr)
        else
            open(newunit = file_unit, file="grid_interface.dat", iostat = ierr)
        end if
        
        if(ierr .ne. 0) then !unable to open file
            res = ierr
            return
        end if
        
        !write things...
        associate( np => this%n_points, &
                    s => (this%start + this%volume_element), &
                    e => (this%start + this%n_points * this%volume_element) )
        
        if(ierr == 0) write(file_unit,"('========== INTERFACE_GRID ==========')", iostat = ierr)
        if(ierr == 0) write(file_unit,"('| d | n_po |   start   |    end    |')", iostat = ierr)
        if(ierr == 0) write(file_unit,"('| X | ',I4,' | ',F9.3,' | ',F9.3,' |')", iostat = ierr) np(1), s(1), e(1)
        if(ierr == 0) write(file_unit,"('| Y | ',I4,' | ',F9.3,' | ',F9.3,' |')", iostat = ierr) np(2), s(2), e(2)
        if(ierr == 0) write(file_unit,"('| Z | ',I4,' | ',F9.3,' | ',F9.3,' |')", iostat = ierr) np(3), s(3), e(3)
        if(ierr == 0) write(file_unit,"('====================================')", iostat = ierr)
        
        end associate

        if(ierr == 0) close(file_unit)
        
        res = ierr
        
    end function write_grid_interface

    function open_bin_file(this, read_only, must_exist, name) result(res)
        class(instantaneous_surface_type), intent(inout) :: this
        character(*), intent(in), optional :: name
        logical, intent(in), optional :: read_only, must_exist
        character(128) :: file_name
        integer :: res
        logical :: is_open
        character(20) :: act, stat
        
        inquire(this%instasurf_bin_unit, opened = is_open)
        if(is_open) close(this%instasurf_bin_unit)
        
        file_name = "interface.bin"
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
            
        open(newunit = this%instasurf_bin_unit, &
                file = trim(adjustl(file_name)), &
                status = stat, &
                action = act, &
                form = "unformatted", &
                access = "stream", &
                convert = "little_endian", &
                iostat = res)
    end function open_bin_file

    function open_xyz_file(this, name) result(res)
        class(instantaneous_surface_type), intent(inout) :: this
        character(*), intent(in), optional :: name
        logical :: is_open
        integer :: res
        character(128) :: file_name
        
        inquire(this%instasurf_xyz_unit, opened = is_open)
        if(is_open) close(this%instasurf_xyz_unit)

        file_name = "interface.xyz"
        
        if(present(name)) then
            file_name = trim(adjustl(name))
        end if
        
        open(newunit = this%instasurf_xyz_unit, &
                file = trim(adjustl(file_name)), &
                iostat = res)
    end function open_xyz_file

    function write_bin_frame(this) result(res)
        class(instantaneous_surface_type), intent(inout) :: this
        logical :: is_open
        integer :: res
        
        inquire(this%instasurf_bin_unit, opened = is_open)
        res = -1; if(.not. is_open) return

        res = 0

        if(res == 0) write(this%instasurf_bin_unit, iostat = res) instantaneous_surface_magic
        if(res == 0) write(this%instasurf_bin_unit, iostat = res) instantaneous_surface_bin_version
        if(res == 0) write(this%instasurf_bin_unit, iostat = res) this%n_points
        if(res == 0) write(this%instasurf_bin_unit, iostat = res) this%volume_element
        if(res == 0) write(this%instasurf_bin_unit, iostat = res) this%start
        if(res == 0) write(this%instasurf_bin_unit, iostat = res) this%up_mesh
        if(res == 0) write(this%instasurf_bin_unit, iostat = res) this%bot_mesh
        
        !results in iostat
    end function write_bin_frame

    function write_xyz_frame(this) result(res)
        class(instantaneous_surface_type), intent(inout) :: this
        logical :: is_open
        integer :: res, i, j
        
        inquire(this%instasurf_xyz_unit, opened = is_open)
        res = -1; if(.not. is_open) return
        
        res = 0
        
        if(res == 0) write(this%instasurf_xyz_unit, "(I0)", iostat = res) (size(this%bot_mesh) + size(this%up_mesh)) !number of points
        if(res == 0) write(this%instasurf_xyz_unit, "(A,I0,A,I0,A)", iostat = res) "bot (", size(this%bot_mesh), ") and up (", size(this%up_mesh), ") interface mesh"

        if(res == 0) then
            do i = 1, size(this%bot_mesh,1)
                do j = 1, size(this%bot_mesh,2)
                    associate( x => this%start(1) + i*this%volume_element(1), &
                                y => this%start(2) + j*this%volume_element(2), &
                                z => this%bot_mesh(i,j) )
                    
                    write(this%instasurf_xyz_unit, "(A, 3F8.3)", iostat = res) "P", x, y, z      
                    if(res .ne. 0) return

                    end associate
                end do
            end do

            do i = 1, size(this%up_mesh,1)
                do j = 1, size(this%up_mesh,2)
                    associate( x => this%start(1) + i*this%volume_element(1), &
                                y => this%start(2) + j*this%volume_element(2), &
                                z => this%up_mesh(i,j) )
                    
                    write(this%instasurf_xyz_unit, "(A, 3F8.3)", iostat = res) "P", x, y, z      
                    if(res .ne. 0) return

                    end associate
                end do
            end do
        end if

    end function write_xyz_frame

    function read_next(this) result(res)
        class(instantaneous_surface_type), intent(inout) :: this
        logical :: is_open
        integer :: res
        integer(int32) :: verify

        inquire(this%instasurf_bin_unit, opened = is_open)
        res = -1; if(.not. is_open) return
        
        res = 0
        
        if(res == 0) then !read and verify magic
            read(this%instasurf_bin_unit, iostat = res) verify
            if(res .ne. 0) return
            if(verify .ne. instantaneous_surface_magic) then
                write(output_unit,*) "ERROR: instantaneous surface magic is not correct"
                res = -2
                return
            end if
        end if

        if(res == 0) then !read and verify version
            read(this%instasurf_bin_unit, iostat = res) verify
            if(res .ne. 0) return
            if(verify .ne. instantaneous_surface_bin_version) then
                write(output_unit,*) "ERROR: instantaneous surface version is not correct"
                res = -3
                return
            end if
        end if

        if(res == 0) then !read n_points
            read(this%instasurf_bin_unit, iostat = res) this%n_points
            if(res .ne. 0) return

            if(allocated(this%up_mesh)) then
                !verify size against n_points
                if( (size(this%up_mesh,1) .ne. this%n_points(1)) .or. &
                    (size(this%up_mesh,1) .ne. this%n_points(1)) ) then
                    !missmatch - deallocate, reallocate
                    deallocate(this%up_mesh)
                    allocate(this%up_mesh(this%n_points(1), this%n_points(2)))
                end if
            else
                !allocate up_mesh
                allocate(this%up_mesh(this%n_points(1), this%n_points(2)))
            end if

            if(allocated(this%bot_mesh)) then
                !verify size against n_points
                if( (size(this%bot_mesh,1) .ne. this%n_points(1)) .or. &
                    (size(this%bot_mesh,2) .ne. this%n_points(2)) ) then
                    !missmatch - deallocate, reallocate
                    deallocate(this%bot_mesh)
                    allocate(this%bot_mesh(this%n_points(1), this%n_points(2)))
                end if
            else
                !allocate bot_mesh
                allocate(this%bot_mesh(this%n_points(1), this%n_points(2)))
            end if
        end if

        if(res == 0) read(this%instasurf_bin_unit, iostat = res) this%volume_element
        if(res == 0) read(this%instasurf_bin_unit, iostat = res) this%start
        if(res == 0) read(this%instasurf_bin_unit, iostat = res) this%up_mesh
        if(res == 0) read(this%instasurf_bin_unit, iostat = res) this%bot_mesh

    end function read_next

    subroutine close_bin_file(this)
        class(instantaneous_surface_type), intent(inout) :: this
        logical :: is_open
        inquire(this%instasurf_bin_unit, opened = is_open)
        if(is_open) close(this%instasurf_bin_unit)
    end subroutine close_bin_file 
    
    subroutine close_xyz_file(this)
        class(instantaneous_surface_type), intent(inout) :: this
        logical :: is_open
        inquire(this%instasurf_xyz_unit, opened = is_open)
        if(is_open) close(this%instasurf_xyz_unit)
    end subroutine close_xyz_file 

end module
