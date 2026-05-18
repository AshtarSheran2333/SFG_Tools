module BOXDATA
    use iso_fortran_env
    use PRETTY_PRINT
    use UTILS_ERROR
    implicit none

#include "utils_error_macros.h"
    
    type boxdata_type
        !PRIVATE variables
        integer, private :: fileUnit

        character(len=3), private :: polarization = "SSP"	
        
        !PUBLIC variables (boxdata parameters) - TODO No time to write getters - it will be public
        real(real64), dimension(3) :: box_dimensions = (/-1,-1,-1/),& !Angstrom
                                      box_corner = (/0.0, 0.0, 0.0/),& !Angstrom
                                      interface_volume_element = (/0.5d0, 0.5d0, 0.25d0/) !Angstrom

        real(real64), dimension(:), allocatable :: layers_limits !Angstrom
                                                            
        real(real64) :: FREQ = 4000,& !cm-1
                        DFREQ = 1,& !cm-1
                        DT = -1,& !fs
                        CORRLEN = 5,& !ps
                        FILTER = 1.5,& !ps
                        TEMPERATURE = 300,& !K
                        DENSITY_R_START = -5.d0,& !min distance from instasurf A
                        DENSITY_R_END = 80.d0,& !max distance from instasurf A
                        LIQUID_BULK_NUMBER_DENSITY = 0.03336,& !expected liquid heavy atom bulk number density -/Angstrom^3
                        COARSE_GRAINING_LENGTH = 2.4 !A
        
        integer(int64) :: NSTEP=-1,& !number of steps -
                          INTERFACE_SKIP=1,& !number of frames skipped when evaluating IS -
                          CROSS_SKIP=1,& !number time beginnings to skip for crossterms -
                          SELF_SKIP=1,& !number of time beginnings to skip for selfterms -
                          INTERFACE_PUSHBACK = 100000,& !-
                          HBHIST_SKIP=1,& !number of frames skipped when evaluating hydrogen bonds -
                          DENSITY_N_POINTS = 1000,&
                          !polarization setup variables
                          P = 1, & !chi2_{pqr}
                          Q = 1, & !chi2_{pqr}
                          R = 3 !chi2_{pqr}

        integer(int32), dimension(:), allocatable :: LIQUID_HEAVY_ATOMS !list of indices of liquid heavy atoms
        
    contains

        procedure, public :: read_boxdata
        procedure, public :: get_maxlag !TODO probably handled elsewhere

        procedure, private :: print_recap
        procedure, private :: read_box_dimensions
        procedure, private :: read_box_corner
        procedure, private :: read_layers_limits
        procedure, private :: read_interface_volume_element
        procedure, private :: read_freq
        procedure, private :: read_dfreq
        procedure, private :: read_dt
        procedure, private :: read_corrlen
        procedure, private :: read_filter
        procedure, private :: read_temperature
        procedure, private :: read_density_r_start
        procedure, private :: read_density_r_end
        procedure, private :: read_density_n_points
        procedure, private :: read_liquid_bulk_number_density
        procedure, private :: read_coarse_graining_length
        procedure, private :: read_nstep
        procedure, private :: read_interface_skip
        procedure, private :: read_cross_skip
        procedure, private :: read_self_skip
        procedure, private :: read_interface_pushback
        procedure, private :: read_hbhist_skip
        procedure, private :: read_liquid_heavy_atoms
        procedure, private :: read_polarization

        procedure, private :: post_read_values_validation
    end type boxdata_type
    
contains

function post_read_values_validation(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res
    real(real64) :: temp_real
    
    res = 0
    
    !check that all the values are nonzero positive
    if(any(this%box_dimensions <= 0)) then
        res = res - 1
        write(error_unit,f_line) "BOXDATA ERROR: $BOX_DIMENSIONS was not set or contains invalid values (<= 0)"
    end if

    !no check
    !this%box_corner

    !check that all the values are nonzero positive
    if(any(this%interface_volume_element <= 0)) then
        res = res - 1
        write(error_unit,f_line) "BOXDATA ERROR: $INTERFACE_VOLUME_ELEMENT contains invalid values (<= 0)"
    end if

    !no check 
    !this%layers_limits
    !TODO check ascending order!!!
                                                            
    !check nonzero, positive 
    if(this%FREQ <= 0.0) then
        write(output_unit,f_line) "BOXDATA WARNING: $FREQ - invalid value, reset to default (4000)"
        this%FREQ = 4000.0
    end if

    !check nonzero, positive 
    if(this%DFREQ <= 0.0) then
        write(output_unit,f_line) "BOXDATA WARNING: $DFREQ - invalid value, reset to default (1)"
        this%DFREQ = 1.0
    end if

    !check nonzero, positive 
    if(this%DT <= 0.0) then
        res = res - 1
        write(error_unit,f_line) "BOXDATA ERROR: $DT was not set or contains invalid value (<= 0)"
    end if

    !check nonzero, positive 
    if(this%CORRLEN <= 0.0) then
        write(output_unit,f_line) "BOXDATA WARNING: $CORRLEN - invalid value, reset to default (5)"
        this%CORRLEN = 5.0
    end if

    !check nonzero, positive 
    if(this%FILTER <= 0.0) then
        write(output_unit,f_line) "BOXDATA WARNING: $FILTER - invalid value, reset to default (1.5)"
        this%FILTER = 1.5
    end if

    !check nonzero, positive 
    if(this%TEMPERATURE <= 0.0) then
        res = res - 1
        write(error_unit,f_line) "BOXDATA ERROR: $TEMPERATURE contains invalid value (<= 0)"
    end if

    !check that start is < end 
    if(this%DENSITY_R_START > this%DENSITY_R_END) then
        write(output_unit,f_line) "BOXDATA WARNING: $DENSITY_R_START > $DENSITY_R_END - swapping values"
        temp_real = this%DENSITY_R_START
        this%DENSITY_R_START = this%DENSITY_R_END
        this%DENSITY_R_END = temp_real
    end if

    !check nonzero, positive 
    if(this%LIQUID_BULK_NUMBER_DENSITY <= 0.0) then
        res = res - 1
        write(error_unit,f_line) "BOXDATA ERROR: $LIQUID_BULK_NUMBER_DENSITY contains invalid value (<= 0)"
    end if

    !check nonzero, positive
    if(this%COARSE_GRAINING_LENGTH <= 0.0) then
        res = res - 1
        write(error_unit,f_line) "BOXDATA ERROR: $COARSE_GRAINING_LENGTH contains invalid value (<= 0)"
    end if

    !check nonzero, positive
    if(this%NSTEP <= 0.0) then
        res = res - 1
        write(error_unit,f_line) "BOXDATA ERROR: $NSTEP was not set or contains invalid value (<= 0)"
    end if

    !check <= 0 ---> 1
    if(this%INTERFACE_SKIP <= 0.0) then
        write(output_unit,f_line) "BOXDATA WARNING: $INTERFACE_SKIP - invalid value, reset to default (1)"
        this%INTERFACE_SKIP = 1
    end if

    !check <= 0 ---> 1
    if(this%CROSS_SKIP <= 0.0) then
        write(output_unit,f_line) "BOXDATA WARNING: $CROSS_SKIP - invalid value, reset to default (1)"
        this%CROSS_SKIP = 1
    end if

    !check <= 0 ---> 1
    if(this%SELF_SKIP <= 0.0) then
        write(output_unit,f_line) "BOXDATA WARNING: $SELF_SKIP - invalid value, reset to default (1)"
        this%SELF_SKIP = 1
    end if

    !check nonzero, positive
    if(this%INTERFACE_PUSHBACK <= 0.0) then
        write(output_unit,f_line) "BOXDATA WARNING: $INTERFACE_PUSHBACK - invalid value, reset to default (100000)"
        this%INTERFACE_PUSHBACK = 100000
    end if
    
    !check <= 0 ---> 1
    if(this%HBHIST_SKIP <= 0.0) then
        write(output_unit,f_line) "BOXDATA WARNING: $HBHIST_SKIP - invalid value, reset to default (1)"
        this%HBHIST_SKIP = 1
    end if

    !check nonzero, positive -> default
    if(this%DENSITY_N_POINTS <= 0.0) then
        write(output_unit,f_line) "BOXDATA WARNING: $DENSITY_N_POINTS - invalid value, reset to default (1000)"
        this%DENSITY_N_POINTS = 1000
    end if

    !no check
    !this%LIQUID_HEAVY_ATOMS

end function post_read_values_validation

function read_box_dimensions(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: i, ierr
    integer :: res

    res = 0;
    
    do i=1,3
        read(this%fileUnit,*, iostat = ierr) this%box_dimensions(i)
        if(ierr .ne. 0) then
            write(error_unit,"(A,I0,A)") "BOXDATA ERROR: $BOX_DIMENSIONS(", i, ") does not contain proper data"
            res = -1
        end if
    end do
end function read_box_dimensions

function read_box_corner(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: i, ierr
    integer :: res

    res = 0;
    
    do i=1,3
        read(this%fileUnit,*, iostat = ierr) this%box_corner(i)
        if(ierr .ne. 0) then
            write(error_unit,"(A,I0,A)") "BOXDATA ERROR: $BOX_CORNER(", i, ") does not contain proper data"
            res = -1
        end if
    end do
end function read_box_corner

function read_layers_limits(this) result(res)
    implicit none
    class(boxdata_type) :: this
    character(len = 256) :: line
    integer :: i, ierr
    integer :: res
    real(real64) :: limit
    real(real64), dimension(:), allocatable :: temp

    if(allocated(this%layers_limits)) deallocate(this%layers_limits)

    res = 0

    do i=1,7
        read(this%fileUnit,*, iostat = ierr) line
        if(ierr .ne. 0) exit !if unable to read, it is the end of the file

        read(line,*, iostat = ierr) limit
        if(ierr .eq. 0) then
            !append the layers_limits
            allocate(temp(size(this%layers_limits)+1))
            temp(1:size(this%layers_limits)) = this%layers_limits
            temp(size(temp)) = limit
            call move_alloc(from=temp, to=this%layers_limits)
        else
            !no more numbers to read
            backspace(this%fileUnit, iostat = ierr) !put back the line
            error_io_check(ierr, "BOXDATA reading fatal error") !should not happen
            exit
        end if
    end do
end function read_layers_limits

function read_interface_volume_element(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: i, ierr
    integer :: res

    res = 0

    do i=1,3
        read(this%fileUnit,*, iostat = ierr) this%interface_volume_element(i)
        if(ierr .ne. 0) then
            write(error_unit,"(a,i1,a)") "BOXDATA ERROR: $INTERFACE_VOLUME_ELEMENT(", i, ") does not contain proper data"
            res = -1
        end if
    end do
end function read_interface_volume_element

function read_freq(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%FREQ
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $FREQ does not contain proper data"
        res = -1
    end if
end function read_freq

function read_dfreq(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*) this%DFREQ
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $DFREQ does not contain proper data"
        res = -1
    end if
end function read_dfreq

function read_dt(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%DT
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $DT does not contain proper data"
        res = -1
    end if
end function read_dt

function read_corrlen(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%CORRLEN
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $CORRLEN does not contain proper data"
        res = -1
    end if
end function read_corrlen

function read_filter(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%FILTER
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $FILTER does not contain proper data"
        res = -1
    end if
end function read_filter

function read_temperature(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%TEMPERATURE
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $TEMPERATURE does not contain proper data"
        res = -1
    end if
end function read_temperature

function read_density_r_start(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0

    read(this%fileUnit,*, iostat = ierr) this%density_r_start
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $DENSITY_R_START does not contain proper data"
        res = -1
    end if
end function read_density_r_start

function read_density_r_end(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%density_r_end
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $DENSITY_R_END does not contain proper data"
        res = -1
    end if
end function read_density_r_end

function read_density_n_points(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0

    read(this%fileUnit,*, iostat = ierr) this%density_n_points
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $DENSITY_N_POINTS does not contain proper data"
        res = -1
    end if
end function read_density_n_points

function read_liquid_bulk_number_density(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%liquid_bulk_number_density
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $LIQUID_BULK_NUMBER_DENSITY does not contain proper data"
        res = -1
    end if
end function read_liquid_bulk_number_density

function read_coarse_graining_length(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%coarse_graining_length
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $COARSE_GRAINING_LENGTH does not contain proper data"
        res = -1
    end if
end function read_coarse_graining_length

function read_nstep(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%NSTEP
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $NSTEP does not contain proper data"
        res = -1
    end if
end function read_nstep

function read_interface_skip(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%INTERFACE_SKIP
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $INTERFACE_SKIP does not contain proper data"
        res = -1
    end if
end function read_interface_skip

function read_cross_skip(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%cross_skip
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $CROSS_SKIP does not contain proper data"
        res = -1
    end if
end function read_cross_skip

function read_self_skip(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%SELF_SKIP
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $SELF_SKIP does not contain proper data"
        res = -1
    end if
end function read_self_skip

function read_interface_pushback(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%INTERFACE_PUSHBACK
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $INTERFACE_PUSHBACK does not contain proper data"
        res = -1
    end if
end function read_interface_pushback

function read_hbhist_skip(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%HBHIST_SKIP
    if(ierr .ne. 0) then
        write(error_unit,f_line) "BOXDATA ERROR: $HBHIST_SKIP does not contain proper data"
        res = -1
    end if
end function read_hbhist_skip

function read_polarization(this) result(res)
    implicit none
    class(boxdata_type) :: this
    character(len=256) :: line
    integer :: res, ierr, i

    res = 0

    read(this%fileUnit,'(A)') line
    line = trim(line)
    line = adjustl(line)
    
    !uppercase
    do i = 1, len(line)
        if(line(i:i) >= 'a' .and. line(i:i) <= 'z') then
            line(i:i) = achar(iachar(line(i:i)) - 32)
        end if
    end do
    
    select case(line)
        case ('SSP')
            this%P=1
            this%Q=1
            this%R=3
            this%polarization = "SSP"	

        case ('PPP')
            this%P=3
            this%Q=3
            this%R=3
            this%polarization = "PPP"	
        
        case default
            write(error_unit,"(A,A,A)") "BOXDATA ERROR: $POLARIZATION ", line, " is not supported."
            res = -1
    end select
end function read_polarization

function read_liquid_heavy_atoms(this) result(res)
    implicit none
    class(boxdata_type) :: this
    character(256) :: line
    integer :: res, ierr, n, sz
    integer(int32), dimension(:), allocatable :: temp
    integer(int32), dimension(32) :: numbers
    
    res = 0
    
    ierr = 0
    do while(ierr == 0)
        read(this%fileUnit, "(A)", iostat = ierr) line
        if(ierr .ne. 0) exit

        numbers = 0
        read(line, *, iostat = ierr) (numbers(n), n = 1, size(numbers))
        if(n == 1) then !read until line does not contain a single number
            backspace(this%fileUnit)
            exit
        end if
        n = n - 1
        
        if(.not. allocated(this%liquid_heavy_atoms)) then
            allocate(this%liquid_heavy_atoms(n))
            this%liquid_heavy_atoms = numbers(1:n)
        else
            sz = n + size(this%liquid_heavy_atoms)
            allocate(temp(sz))
            sz = size(this%liquid_heavy_atoms)
            temp(1:sz) = this%liquid_heavy_atoms
            temp(sz+1:) = numbers(1:n)
            call move_alloc(from = temp, to = this%liquid_heavy_atoms)
        end if
        ierr = 0
    end do

    res = 0 !res is allways 0...
end function read_liquid_heavy_atoms

subroutine read_boxdata(this)
    class(boxdata_type) :: this
    character(len = 256) :: line
    integer :: ierr,&
               i,&
               error_count

    error_count = 0
    
    write(output_unit,f_line) heading(flat_pattern, "Opening BOXDATA")
    write(output_unit,f_line) ""
    
    open(newunit = this%fileUnit, file = "BOXDATA", status = 'old', iostat = ierr)
    error_io_check(ierr, "Problem opening BOXDATA file")
    
    read(this%fileUnit,'(A)', iostat = ierr) line
    
    !evaluate all the lines in BOXDATA
    do while (ierr .eq. 0)
            
        line = trim(adjustl(line))
        
        select case(line)
        
            case("$BOX_DIMENSIONS") !box dimensions
                error_count = error_count + this%read_box_dimensions()
                
            case("$BOX_CORNER") !box dimensions
                error_count = error_count + this%read_box_corner()

            case("$INTERFACE_VOLUME_ELEMENT") !box dimensions
                error_count = error_count + this%read_interface_volume_element()
                
            case("$LAYERS_LIMITS") !Z coordinates of layers
                error_count = error_count + this%read_layers_limits()
                
            case("$FREQ") !frequency range for the fourier transform
                error_count = error_count + this%read_freq()
                
            case("$DFREQ") !frequency range for the fourier transform
                error_count = error_count + this%read_dfreq()

            case("$DT") !timestep in fs
                error_count = error_count + this%read_dt()

            case("$FILTER") !number of skipped frames while calculating correlation function
                error_count = error_count + this%read_filter()
            
            case("$CORRLEN") !length of correlation function in ps
                error_count = error_count + this%read_corrlen()
               
            case("$TEMPERATURE")
                error_count = error_count + this%read_temperature()

            case("$DENSITY_R_START")
                error_count = error_count + this%read_density_r_start()

            case("$DENSITY_R_END")
                error_count = error_count + this%read_density_r_end()

            case("$LIQUID_BULK_NUMBER_DENSITY")
                error_count = error_count + this%read_liquid_bulk_number_density()

            case("$COARSE_GRAINING_LENGTH")
                error_count = error_count + this%read_coarse_graining_length()
                
            case("$NSTEP") !length of trajectory
                error_count = error_count + this%read_nstep()
                
            case("$INTERFACE_SKIP") !number of skipped frames to calculate interface
                error_count = error_count + this%read_interface_skip()
                
            case("$CROSS_SKIP") !number of skipped frames while calculating correlation function
                error_count = error_count + this%read_cross_skip()
                    
            case("$SELF_SKIP") !number of skipped frames while calculating correlation function
                error_count = error_count + this%read_self_skip()
                    
            case("$INTERFACE_PUSHBACK") !number of skipped frames while calculating correlation function
                error_count = error_count + this%read_interface_pushback()

            case("$HBHIST_SKIP") !number of skipped frames while calculating correlation function
                error_count = error_count + this%read_hbhist_skip()

            case("$DENSITY_N_POINTS") !number of skipped frames while calculating correlation function
                error_count = error_count + this%read_density_n_points()

            case("$POLARIZATION")
                error_count = error_count + this%read_polarization()

            case("$LIQUID_HEAVY_ATOMS")
                error_count = error_count + this%read_liquid_heavy_atoms() !list of heavy atoms in the liquid
                
            case("") !empty line
                
            case default
                if(index(line, '#') .ne. 1) then
                    write(output_unit,"(A,' ',A)") "BOXDATA WARNING: skipping invalid line:", trim(line)
                end if
            end select
                
        read(this%fileUnit,'(A)', IOSTAT = ierr) line
        
    end do

    close(this%fileUnit)
    
    error_count = error_count + this%post_read_values_validation()
    
    call this%print_recap()

    if(error_count < 0) then !if error in reading the file
        error_stop("BOXDATA file contains some errors, fix them please")
    else
        write(output_unit,f_line) heading(flat_pattern, "BOXDATA - DONE")
        write(output_unit,f_line) ""
    end if
    
end subroutine read_boxdata

!!!returns number of samples of correlation fucntion based on DT and CORRLEN (ps)
integer*8 function get_maxlag(this)
    class(boxdata_type) :: this
    
    if(this%corrlen == 0) then
        get_maxlag = this%nstep - 1
    else if(nint(1000*this%CORRLEN/this%dt) >= this%nstep) then
        get_maxlag = this%nstep - 1
    else
        get_maxlag = nint(1000*this%CORRLEN/this%DT)
    end if

end function get_maxlag

subroutine print_recap(this)
    class(boxdata_type) :: this
    integer :: i
    
    !TODO update this function when everything is settled
    write(output_unit,f_line) ""
    write(output_unit,f_line) "BOXDATA RECAP:"
    write(output_unit,f_line) ""

    !arrays
    write(output_unit,"(A40,5X,F10.3)") "BOX_DIMENSIONS", this%box_dimensions(1)
    write(output_unit,"(A40,5X,F10.3)") "", this%box_dimensions(2)
    write(output_unit,"(A40,5X,F10.3)") "", this%box_dimensions(3)
    write(output_unit,"(A40,5X,F10.3)") "BOX_CORNER", this%box_corner(1)
    write(output_unit,"(A40,5X,F10.3)") "", this%box_corner(2)
    write(output_unit,"(A40,5X,F10.3)") "", this%box_corner(3)
    write(output_unit,"(A40,5X,F10.3)") "INTERFACE_VOLUME_ELEMENT", this%interface_volume_element(1)
    write(output_unit,"(A40,5X,F10.3)") "", this%interface_volume_element(2)
    write(output_unit,"(A40,5X,F10.3)") "", this%interface_volume_element(3)
    write(output_unit,"(A40,5X,F10.3)") "LAYERS_LIMITS", this%layers_limits(1)
    do i=2, size(this%layers_limits)
        write(output_unit,"(A40,5X,F10.3)") "", this%layers_limits(i)
    end do

    !reals
    write(output_unit,"(A40,5X,F10.3)") "FREQ", this%FREQ
    write(output_unit,"(A40,5X,F10.3)") "DFREQ", this%DFREQ
    write(output_unit,"(A40,5X,F10.3)") "DT", this%DT
    write(output_unit,"(A40,5X,F10.3)") "CORRLEN", this%CORRLEN
    write(output_unit,"(A40,5X,F10.3)") "FILTER", this%FILTER
    write(output_unit,"(A40,5X,F10.3)") "TEMPERATURE", this%TEMPERATURE
    write(output_unit,"(A40,5X,F10.3)") "DENSITY_R_START", this%DENSITY_R_START
    write(output_unit,"(A40,5X,F10.3)") "DENSITY_R_END", this%DENSITY_R_END
    write(output_unit,"(A40,5X,F10.3)") "LIQUID_BULK_NUMBER_DENSITY", this%LIQUID_BULK_NUMBER_DENSITY
    write(output_unit,"(A40,5X,F10.3)") "COARSE_GRAINING_LENGTH", this%COARSE_GRAINING_LENGTH

    !integers
    write(output_unit,"(A40,5X,I10)") "NSTEP", this%NSTEP
    write(output_unit,"(A40,5X,I10)") "INTERFACE_SKIP", this%INTERFACE_SKIP
    write(output_unit,"(A40,5X,I10)") "CROSS_SKIP", this%CROSS_SKIP
    write(output_unit,"(A40,5X,I10)") "SELF_SKIP", this%SELF_SKIP
    write(output_unit,"(A40,5X,I10)") "INTERFACE_PUSHBACK", this%INTERFACE_PUSHBACK
    write(output_unit,"(A40,5X,I10)") "HBHIST_SKIP", this%HBHIST_SKIP
    write(output_unit,"(A40,5X,I10)") "DENSITY_N_POINTS", this%DENSITY_N_POINTS

    !chars
    write(output_unit,"(A40,5X,A10)") "POLARIZATION", this%polarization

    !logicals

    write(output_unit,f_line) ""

end subroutine print_recap

end module BOXDATA
    
