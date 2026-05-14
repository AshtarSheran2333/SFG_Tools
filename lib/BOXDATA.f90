module BOXDATA
    use iso_fortran_env
    use UTILS_ERROR
    implicit none

#include "utils_error_macros.h"
    
    type boxdata_type
        !PRIVATE variables
        integer, private :: fileUnit,&
                            layersCount = 1

        character(len=3), private :: polarization = "SSP"	
        
        !PUBLIC variables (boxdata parameters) - TODO No time to write getters - it will be public
        real(real64), dimension(3) :: box_dimensions = (/-1,-1,-1/),& !Angstrom
                                      box_corner = (/0.0, 0.0, 0.0/),& !Angstrom
                                      interface_volume_element = (/0.5d0, 0.5d0, 0.25d0/) !Angstrom
        
        real(real64), dimension(7) :: layers_limits = (/0,0,0,0,0,0,0/) !Angstrom
                                                            
        
        real(real64) :: FREQ = 4000,& !cm-1
                        DFREQ = 1,& !cm-1
                        DT = -1,& !fs
                        CORRLEN = 0,& !ps
                        FILTER = 0,& !ps
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
        procedure, public :: get_layersCount

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
    end type boxdata_type
    
contains

function read_box_dimensions(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: i, ierr
    integer :: res

    res = 0;
    
    do i=1,3
        read(this%fileUnit,*, iostat = ierr) this%box_dimensions(i)
        if(ierr .ne. 0) then
            print"(A,I0,A)", "BOXDATA ERROR: $BOX_DIMENSIONS(", i, ") does not contain proper data"
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
            print"(A,I0,A)", "BOXDATA ERROR: $BOX_CORNER(", i, ") does not contain proper data"
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

    !TODO make it better

    res = 0

    do i=1,7
        read(this%fileUnit,*, iostat = ierr) line
        if(ierr .ne. 0) then
            print*, "problem when reading BOXDATA $LAYERS_LIMITS is probably the last line of the file"
        end if
        read(line,*, iostat = ierr) this%layers_limits(i)
        if(ierr .eq. 0) then
            this%layersCount = this%layersCount + 1
        else
            !no more numbers to read
            backspace(this%fileUnit, iostat = ierr)
            if(ierr .ne. 0) then
                print*, "backspace problem when reading $LAYERS_LIMITS"
                stop
            end if
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
            print"(a,i1,a)", "BOXDATA ERROR: $INTERFACE_VOLUME_ELEMENT(", i, ") does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $FREQ does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $DFREQ does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $DT does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $CORRLEN does not contain proper data"
        res = -1
    end if
end function read_corrlen

function read_filter(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    if(this%filter == 0) then
        read(this%fileUnit,*, iostat = ierr) this%FILTER
        if(ierr .ne. 0) then
            print"(a)", "BOXDATA WARNING: $FILTER does not contain proper data"
            print*, "Filter will be set to default value"
            print*, ""
        end if
    else
        print*, "REMINDER:"
        print*, "$FILTER was modified by the program"
        print*, "$FILTER value is set to:", this%filter
    end if
end function read_filter

function read_temperature(this) result(res)
    implicit none
    class(boxdata_type) :: this
    integer :: res, ierr

    res = 0
    
    read(this%fileUnit,*, iostat = ierr) this%TEMPERATURE
    if(ierr .ne. 0) then
        print"(a)", "BOXDATA ERROR: $TEMPERATURE does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $DENSITY_R_START does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $DENSITY_R_end does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $DENSITY_N_POINTS does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $INTERFACE_DENSITY does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $COARSE_GRAINING_LENGTH does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $NSTEP does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $INTERFACE_SKIP does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $CROSS_SKIP does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $SELF_SKIP does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $INTERFACE_PUSHBACK does not contain proper data"
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
        print"(a)", "BOXDATA ERROR: $HBHIST_SKIP does not contain proper data"
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
            print*, "BOXDATA ERROR: polarization ", line, " is not supported."
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
        if(n == 1) then
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

    res = 0
end function read_liquid_heavy_atoms

subroutine read_boxdata(this)
    class(boxdata_type) :: this
    character(len = 256) :: line
    integer :: ierr,&
               i,&
               error_count
    
    open(newunit = this%fileUnit, file = "BOXDATA", status = 'old', iostat = ierr)
    print*, "Opening BOXDATA..."
    print*, ""
    
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
                print*, "BOXDATA REMAINDER: skipping ", trim(line), " - not a valid line"
                
            end select
                
        read(this%fileUnit,'(A)', IOSTAT = ierr) line
        
    end do
    
    !TODO post read check
    
    if(error_count < 0) then !if error in reading the file
        error_stop("BOXDATA file contains some errors, fix them please")
    else
        print*, "reading of the BOXDATA file - OK"
        print*, ""
    end if
    
    !checking the input... catching errors in the initial values
    
    do i=1,3
        if(this%box_dimensions(i) <= 0) then
            print"(a,i1,a)", "ERROR: $BOX_DIMENSIONS(", i, ") is NOT positive nonzero value"
!error = .true.
        end if
    end do
    
    do i=1,3
        if(this%layers_limits(i) == 0) then
            print"(a,i1,a)", "WARNING: $LAYERS_LIMITS(", i, ") is zero"
            print*, "This will prevent Binder program to start"
        end if
    end do
    
    do i=1,3
        if(this%interface_volume_element(i) <= 0) then
            print"(a,i1,a)", "ERROR: $INTERFACE_VOLUME_ELEMENT(", i, ") is NOT positive nonzero value"
!error = .true.
        end if
    end do    
    
    if(this%DFREQ >= this%FREQ) then !dfreq is too high
        print*, "ERROR: $DFREQ is higher or equal to $FREQ"
!error = .true.
    end if
    
    if(this%DT <= 0) then
        print*, "ERROR: $DT is not positive nonzero value"
!error = .true.
    end if
    
    if(this%corrlen == 0) then
        print("(A)"), "WARNING:"
        print*, "$CORRLEN was not declared in BOXDATA file"
        print*, "Length of correlation function is set to $NSTEP in BOXDATA"
        print*, "This might lead to long calculation time"
        print*, ""
    else if(int(1000*this%CORRLEN/this%dt) > this%nstep) then
        print("(A)"), "WARNING:"
        print*, "$CORRLEN was longer than length of the simulation"
        print*, "$CORRLEN is set to the length of the simulation (ps): ", this%dt*this%nstep/1000
        print*, "This might lead to long calculation time"
        print*, ""
    end if
    
    if(this%TEMPERATURE <= 0) then
        print*, "ERROR: $TEMPERATURE is not positive nonzero value"
!error = .true.
    end if
    
    if(this%density_r_end <= this%density_r_start) then
        print*, "ERROR: $DIPOLE_R_END is lower or equal to $DIPOLE_R_START"
!error = .true.
    end if
    
    if(this%NSTEP <= 0) then
        print("(A)"), "ERROR:"
        print*, "$NSTEP is not positive nonzero value"
!error = .true.
    end if    
    
    if(this%INTERFACE_SKIP < 0) then
        print("(A)"), "WARNING:"
        print*, "$INTERFACE_SKIP is negative number - absolute value will be used"
        this%INTERFACE_SKIP = -this%INTERFACE_SKIP
    end if  

    if(this%INTERFACE_SKIP == 0) then
        this%INTERFACE_SKIP = 1
    end if  
    
    if(this%CROSS_SKIP < 0) then
        print("(A)"), "WARNING:"
        print*, "$CROSS_SKIP is negative number - absolute value will be used"
        this%CROSS_SKIP = -this%CROSS_SKIP
    end if  

    if(this%CROSS_SKIP == 0) then
        this%CROSS_SKIP = 1
    end if  
    
    if(this%SELF_SKIP < 0) then
        print("(A)"), "WARNING:"
        print*, "$SELF_SKIP is negative number - absolute value will be used"
        this%SELF_SKIP = -this%SELF_SKIP
    end if  

    if(this%SELF_SKIP == 0) then
        this%SELF_SKIP = 1
    end if  

    if(this%HBHIST_SKIP < 0) then
        print("(A)"), "WARNING:"
        print*, "$HBHIST_SKIP is negative number - absolute value will be used"
        this%HBHIST_SKIP = -this%HBHIST_SKIP
    end if  

    if(this%HBHIST_SKIP == 0) then
        this%HBHIST_SKIP = 1
    end if  

    if (this%INTERFACE_PUSHBACK <= 0) then
        print("(A)"), "WARNING:"
        print*, "$INTERFACE_PUSHBACK in BOXDATA file is set to 0 or negative number"
        print*, "setting $INTERFACE_PUSHBACK to max value"
        print*, "this may lead to longer calculation time"
        this%INTERFACE_PUSHBACK = 100000
    end if    
    
    if(this%filter <= 0) then
        !filter <= 0 makes no sense -> Simones default version...
        !default filter value by Simone after the change of filter
        print("(A)"), "WARNING:"
        print*, "filter was set to a negative or zero value: ", this%filter
        this%filter = 1.06
        print*, "the filter value was modified to value: ", this%filter, " ps"
    else
        !filter is set according to the value in boxdata (can be overriden by -filter switch in program switches evaluation) 
    end if

    if(error_count < 0) then
        print*, ""
        print*, "BOXDATA file has some ERRORS in the input values, please fix them"
        stop
    else
        print*, "checking the BOXDATA values - OK"
        print*, ""
    end if
    
    print*, "BOXDATA reading finished"
    print*, ""
    
    call this%print_recap()
    print*, ""
    
    close(this%fileUnit)
    
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

integer*8 function get_layersCount(this)
    class(boxdata_type) :: this
    
    get_layersCount = this%layersCount
end function get_layersCount

subroutine print_recap(this)
    class(boxdata_type) :: this
    integer :: i
    
    !TODO update this function when everything is settled
    print*, "BOXDATA RECAP:"
    print*, ""

    !arrays
    print "(A40,5X,F10.3)", adjustl("BOX_DIMENSIONS"), this%box_dimensions(1)
    print "(A40,5X,F10.3)", adjustl(""), this%box_dimensions(2)
    print "(A40,5X,F10.3)", adjustl(""), this%box_dimensions(3)
    print "(A40,5X,F10.3)", adjustl("INTERFACE_VOLUME_ELEMENT"), this%interface_volume_element(1)
    print "(A40,5X,F10.3)", adjustl(""), this%interface_volume_element(2)
    print "(A40,5X,F10.3)", adjustl(""), this%interface_volume_element(3)
    print "(A40,5X,F10.3)", adjustl("LAYERS_LIMITS"), this%layers_limits(1)
    do i = 2, this%layersCount-1
        print "(A40,5X,F10.3)", adjustl(""), this%layers_limits(i)
    end do

    !reals
    print "(A40,5X,F10.3)", adjustl("FREQ"), this%FREQ
    print "(A40,5X,F10.3)", adjustl("DFREQ"), this%DFREQ
    print "(A40,5X,F10.3)", adjustl("DT"), this%DT
    print "(A40,5X,F10.3)", adjustl("CORRLEN"), this%CORRLEN
    print "(A40,5X,F10.3)", adjustl("FILTER"), this%FILTER
    print "(A40,5X,F10.3)", adjustl("TEMPERATURE"), this%TEMPERATURE

    !integers
    print "(A40,5X,I10)", adjustl("NSTEP"), this%NSTEP
    print "(A40,5X,I10)", adjustl("INTERFACE_SKIP"), this%INTERFACE_SKIP
    print "(A40,5X,I10)", adjustl("CROSS_SKIP"), this%CROSS_SKIP
    print "(A40,5X,I10)", adjustl("SELF_SKIP"), this%SELF_SKIP
    print "(A40,5X,I10)", adjustl("INTERFACE_PUSHBACK"), this%INTERFACE_PUSHBACK
    print "(A40,5X,I10)", adjustl("HBHIST_SKIP"), this%HBHIST_SKIP

    !chars
    print "(A40,5X,A10)", adjustl("POLARIZATION"), adjustr(this%polarization)

    !logicals
    

end subroutine print_recap

end module BOXDATA
    
