module INSTANTANEOUS_SURFACE
    use iso_fortran_env
    use FRAME_READERS, only: fr_frame
    use SFG_UTILS, only: pbc_minimum_image, pi
    implicit none
    
    type :: instantaneous_surface_structure
        integer(int32), dimension(3) :: n_points !number of iterations through space
        real(real64), dimension(3) :: volume_element
        real(real64), dimension(3) :: start !corner of the mesh
        real(real64), dimension(:,:), allocatable :: up_mesh
        real(real64), dimension(:,:), allocatable :: bot_mesh
        integer(int32), dimension(:,:), allocatable :: up_index !index of iteration to get the point
        integer(int32), dimension(:,:), allocatable :: bot_index !index of iteration to get the point
    end type instantaneous_surface_structure
    
    integer, parameter :: INSTASURF_BIN = 1, INSTASURF_XYZ = 2
    
    integer, private :: instasurf_bin_unit, instasurf_xyz_unit
    type(instantaneous_surface_structure), allocatable, protected :: instasurf
    
contains

    subroutine instasurf_init(box, volume_element, corner) !TODO work with BOXDATA
        implicit none
        real(real64), dimension(3), intent(in) :: box, volume_element, corner
        
        !wipe
        if(allocated(instasurf)) deallocate(instasurf)
        allocate(instasurf)
        
        !get mesh
        instasurf%n_points = NINT(box/volume_element)

        instasurf%volume_element = box/instasurf%n_points

        instasurf%start = corner - instasurf%volume_element
        
        !allocate containers
        
        associate( i => instasurf%n_points(1), &
                   j => instasurf%n_points(1))

        allocate(instasurf%up_mesh(i,j))
        allocate(instasurf%bot_mesh(i,j))
        allocate(instasurf%up_index(i,j))
        allocate(instasurf%bot_index(i,j))
        
        end associate
    end subroutine instasurf_init

    function instasurf_calculate(atom_selection, graining_len, density_threshold, pushback) result(res)
        implicit none
        !TODO - ugly implementatation - calling the same block twice
        integer(int32), dimension(:), intent(IN) :: atom_selection
        real(real64), intent(IN) :: graining_len, density_threshold
        integer(int32), intent(IN) :: pushback !TODO BOXDATA
        integer(int32) :: res
        
        integer(int32) :: i,j,k,m
        logical :: found_up_interface
        logical :: found_bot_interface
        real(real64), dimension(3) :: pos, prev_pos, diff
        real(real64) :: rho, rhodiff, prev_rhodiff, r

        !$OMP PARALLEL DO &
        !$OMP DEFAULT(NONE) &
        !$OMP SHARED(instasurf, fr_frame, atom_selection, pushback, graining_len, density_threshold) &
        !$OMP PRIVATE(i, j, k, m, found_up_interface, found_bot_interface, pos, prev_pos, diff, rho, rhodiff, prev_rhodiff, r) &
        !$OMP REDUCTION(+:res)
        do i=1,instasurf%n_points(1)
            res = 0
            do j=1,instasurf%n_points(2)
            
                found_up_interface = .false.
                found_bot_interface = .false.
                
                pos = instasurf%start
                rhodiff = density_threshold
    
                !probing from bottom
                do k=max(1,instasurf%bot_index(i,j)-pushback), instasurf%n_points(3)
                    prev_pos = pos
                    pos = instasurf%volume_element * (/i,j,k/) + instasurf%start !TODO volume element from BOXDATA
            
                    !evaluate rho
                    rho = 0
            
                    do m=1, size(atom_selection)
                        associate( atom_pos => fr_frame%positions(:,atom_selection(m)) )
                        diff = atom_pos - pos
                        end associate
                        !pbc correction              
                        diff = pbc_minimum_image(diff, (/10.0_real64,10.0_real64,10.0_real64/)) !TODO BOXDATA CORNER!!! 
                        r = norm2(diff)
                        !cutoff after 3 sigma, the value would be too small, save some calculation time
                        if( r <= 3*graining_len ) then
                            rho = rho + exp(-r**2/(2*graining_len**2))/((2*pi*graining_len**2)**1.5)
                        end if  
                    end do
            
                    !shaping the rho
                    prev_rhodiff = rhodiff
                    rhodiff = abs(density_threshold - rho)
            
                    !if the positive derivative is found -> we have found the interface
                    if( (prev_rhodiff < density_threshold) .and. (rhodiff > prev_rhodiff) ) then
                        found_bot_interface = .true.
                        instasurf%bot_mesh(i,j) = prev_pos(3)
                        instasurf%bot_index(i,j) = k - 1
                        exit
                    end if
                end do
                
                pos = instasurf%start
                rhodiff = density_threshold
    
                !poking rod from top
                do k=min(instasurf%n_points(3), instasurf%up_index(i,j)+pushback),1,-1
                    prev_pos = pos
                    pos = instasurf%volume_element * (/i,j,k/) + instasurf%start
            
                    !evaluate rho
                    rho = 0
            
                    do m=1, size(atom_selection)
                        associate( atom_pos => fr_frame%positions(:,atom_selection(m)) )
                        diff = atom_pos - pos
                        end associate
                        !pbc correction              
                        diff = pbc_minimum_image(diff, (/10.0_real64,10.0_real64,10.0_real64/)) !TODO BOXDATA CORNER!!! 
                        r = norm2(diff)
                        !cutoff after 3 sigma, the value would be too small, save some calculation time
                        if( r <= 3*graining_len ) then
                            rho = rho + exp(-r**2/(2*graining_len**2))/((2*pi*graining_len**2)**1.5)
                        end if  
                    end do

                    !shaping the rho
                    prev_rhodiff = rhodiff
                    rhodiff = abs(density_threshold - rho)
            
                    !if the positive derivative is found -> we have found the interface
                    if( (prev_rhodiff < density_threshold) .and. (rhodiff > prev_rhodiff) ) then
                        found_up_interface = .true.
                        instasurf%up_mesh(i,j) = prev_pos(3)
                        instasurf%up_index(i,j) = k - 1
                        exit
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
    end function instasurf_calculate

    function instasurf_get_distances(point, box, corner) result(distances) !TODO BOXDATA
        implicit none
        real(real64), dimension(3), intent(in) :: point, box, corner
        real(real64), dimension(3) :: diff
        real(real64) :: xgrid, ygrid, xt, yu
        real(real64), dimension(2) :: distances
        integer(int32) :: x, x1, y, y1
        
        diff = pbc_minimum_image(point, box, corner)
        
        !find where the point belongs on the grid
        !X
        x = int((diff(1) - (instasurf%start(1) + instasurf%volume_element(1))) / instasurf%volume_element(1)) + 1
        x1 = x + 1
        !Y
        y = int((diff(2) - (instasurf%start(2) + instasurf%volume_element(2))) / instasurf%volume_element(2)) + 1
        y1 = y + 1

        !deal with the edge cases
        if(x > instasurf%n_points(1)) x = 1
        if(x1 > instasurf%n_points(1)) x1 = 1
        if(y > instasurf%n_points(2)) y = 1
        if(y1 > instasurf%n_points(2)) y1 = 1
        
        xgrid = instasurf%start(1) + x*instasurf%volume_element(1)
        ygrid = instasurf%start(2) + y*instasurf%volume_element(2)
        
        xt = (diff(1)-xgrid)/(instasurf%volume_element(1))
        yu = (diff(2)-ygrid)/(instasurf%volume_element(2))
        
        !interface z
        distances(2) = (1-xt)*(1-yu)*instasurf%up_mesh(x,y) &
                    + xt*(1-yu)*instasurf%up_mesh(x1,y) &
                    + (1-xt)*yu*instasurf%up_mesh(x,y1) &
                    + xt*yu*instasurf%up_mesh(x1,y1)

        distances(1) = (1-xt)*(1-yu)*instasurf%bot_mesh(x,y) &
                    + xt*(1-yu)*instasurf%bot_mesh(x1,y) &
                    + (1-xt)*yu*instasurf%bot_mesh(x,y1) &
                    + xt*yu*instasurf%bot_mesh(x1,y1)
        !distance from the interface
        distances(1) = diff(3) - distances(1)
        distances(2) = distances(2) - diff(3)
    end function

    function instasurf_write_grid_interface(filename) result(res)
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
        associate( np => instasurf%n_points, &
                    s => (instasurf%start + instasurf%volume_element), &
                    e => (instasurf%start + instasurf%n_points * instasurf%volume_element) )
        
        if(ierr == 0) write(file_unit,"('========== INTERFACE_GRID ==========')", iostat = ierr)
        if(ierr == 0) write(file_unit,"('| d | n_po |   start   |    end    |')", iostat = ierr)
        if(ierr == 0) write(file_unit,"('| X | ',I4,' | ',F9.3,' | ',F9.3,' |')", iostat = ierr) np(1), s(1), e(1)
        if(ierr == 0) write(file_unit,"('| Y | ',I4,' | ',F9.3,' | ',F9.3,' |')", iostat = ierr) np(2), s(2), e(2)
        if(ierr == 0) write(file_unit,"('| Z | ',I4,' | ',F9.3,' | ',F9.3,' |')", iostat = ierr) np(3), s(3), e(3)
        if(ierr == 0) write(file_unit,"('====================================')", iostat = ierr)
        
        end associate

        if(ierr == 0) close(file_unit)
        
        res = ierr
        
    end function instasurf_write_grid_interface

    function instasurf_open_bin_file(read_only, must_exist, name) result(res)
        character(*), intent(in), optional :: name
        logical, intent(in), optional :: read_only, must_exist
        character(128) :: file_name
        integer :: res
        logical :: is_open
        character(20) :: act, stat
        
        inquire(instasurf_bin_unit, opened = is_open)
        if(is_open) close(instasurf_bin_unit)
        
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
            
        open(newunit = instasurf_bin_unit, file = trim(adjustl(file_name)), status = stat, action = act, iostat = res)
    end function

    function instasurf_open_xyz_file(name) result(res)
        character(*), intent(in), optional :: name
        logical :: is_open
        integer :: res
        character(128) :: file_name
        
        inquire(instasurf_xyz_unit, opened = is_open)
        if(is_open) close(instasurf_bin_unit)

        file_name = "interface.xyz"
        
        if(present(name)) then
            file_name = trim(adjustl(name))
        end if
        
        open(newunit = instasurf_xyz_unit, file = trim(adjustl(file_name)), iostat = res)
    end subroutine

    subroutine instasurf_write_bin_file()!(name)
        !open the instantaneous surfcace file
        stop "NOT IMPLEMENTED"
    end subroutine

    subroutine instasurf_write_xyz_file()!(name)
        !open the instantaneous surfcace file
        stop "NOT IMPLEMENTED"
    end subroutine

    subroutine instasurf_read_next()!(filetype)
        !read the instantaneous surface from a file
        stop "NOT IMPLEMENTED"
    end subroutine instasurf_read_next

    subroutine instasurf_close_bin_file()!
        !close the instantaneous surface file
        stop "NOT IMPLEMENTED"
    end subroutine instasurf_close_bin_file 
    
    subroutine instasurf_close_xyz_file()!
        !close the instantaneous surface file
        stop "NOT IMPLEMENTED"
    end subroutine instasurf_close_xyz_file 

end module