module DXDRZ_DB
    use iso_fortran_env
    implicit none
    
    type dAdrz_type
        real(real64), dimension(3,3) :: elements
    end type
    
    type dMdrz_type
        real(real64), dimension(3) :: elements
    end type
    
    type dXdrz_db_type
        
        !members kept public to avoid copies... (do not write into them)
        type(dAdrz_type), dimension(:), allocatable :: dAdrz_record
        type(dMdrz_type), dimension(:), allocatable :: dMdrz_record
        integer(int32) :: count
        
    contains
    
        procedure, public :: init
        procedure, public :: read_from_file
    end type

    contains

        subroutine init(this)
            class(dXdrz_db_type), intent(inout) :: this
            
            if(allocated(this%dAdrz_record)) deallocate(this%dAdrz_record)
            allocate(this%dAdrz_record(2))
            if(allocated(this%dMdrz_record)) deallocate(this%dMdrz_record)
            allocate(this%dMdrz_record(2))

            this%count = 2

            !------------------------------------------------------------------------------!
            !                                                                              !
            !   values that does not correspond to anything, isotropic...                  !
            !   where x,y,z are coordinates in the oh_ref                                  !
            !                                                                              !
            !------------------------------------------------------------------------------!
            ! dMdrz(1)   =  0.00 dM(x)/drz
            ! dMdrz(2)   =  0.00 dM(y)/drz
            ! dMdrz(3)   =  1.00 dM(z)/drz
            !
            ! dAdrz(1,1) =  0.50 dA(x,x)/drz
            ! dAdrz(2,2) =  0.50 dA(y,y)/drz
            ! dAdrz(3,3) =  1.00 dA(z,z)/drz
            !
            ! dAdrz(1,2) =  0.00 dA(x,y)/drz
            ! dAdrz(2,1) =  0.00 dA(y,x)/drz
            !
            ! dAdrz(1,3) =  0.00 dA(x,z)/drz
            ! dAdrz(3,1) =  0.00 dA(z,x)/drz
            !
            ! dAdrz(2,3) =  0.00 dA(y,z)/drz
            ! dAdrz(3,2) =  0.00 dA(z,y)/drz

            this%dMdrz_record(2)%elements(:) =   (/  0.00_real64,  0.00_real64,  1.00_real64 /)
            this%dAdrz_record(2)%elements(1,:) = (/  0.50_real64,  0.00_real64,  0.00_real64 /)
            this%dAdrz_record(2)%elements(2,:) = (/  0.00_real64,  0.50_real64,  0.00_real64 /)
            this%dAdrz_record(2)%elements(3,:) = (/  0.00_real64,  0.00_real64,  1.00_real64 /)

            !------------------------------------------------------------------------------!
            !                                                                              !
            !   values from Remi Khatib article for : dM(x,y,z)/dRz and dA(x,y,z)/dRz      !
            !   for water https://doi.org/10.1038/srep24287                                !
            !   where x,y,z are coordinates in the oh_ref                                  !
            !                                                                              !
            !------------------------------------------------------------------------------!
            ! dMdrz(1)   = -0.15 dM(x)/drz
            ! dMdrz(2)   =  0.00 dM(y)/drz
            ! dMdrz(3)   =  2.10 dM(z)/drz
            !
            ! dAdrz(1,1) =  0.40 dA(x,x)/drz
            ! dAdrz(2,2) =  0.53 dA(y,y)/drz
            ! dAdrz(3,3) =  1.56 dA(z,z)/drz
            !
            ! dAdrz(1,2) =  0.00 dA(x,y)/drz
            ! dAdrz(2,1) =  0.00 dA(y,x)/drz
            !
            ! dAdrz(1,3) =  0.02 dA(x,z)/drz
            ! dAdrz(3,1) =  0.02 dA(z,x)/drz
            !
            ! dAdrz(2,3) =  0.00 dA(y,z)/drz
            ! dAdrz(3,2) =  0.00 dA(z,y)/drz

            this%dMdrz_record(2)%elements(:) =   (/ -0.15_real64,  0.00_real64,  2.10_real64 /)
            this%dAdrz_record(2)%elements(1,:) = (/  0.40_real64,  0.00_real64,  0.02_real64 /)
            this%dAdrz_record(2)%elements(2,:) = (/  0.00_real64,  0.53_real64,  0.00_real64 /)
            this%dAdrz_record(2)%elements(3,:) = (/  0.02_real64,  0.00_real64,  1.56_real64 /)
        
        end subroutine init
        
        ! the file is expected to look like this:
        ! #COMMENT
        !
        ! $PARAMETERS
        ! dMdrz(1) dMdrz(2) dMdrz(3)
        ! dAdrz(1,1) dAdrz(1,2) dAdrz(1,3)
        ! dAdrz(2,1) dAdrz(2,2) dAdrz(2,3)
        ! dAdrz(3,1) dAdrz(3,2) dAdrz(3,3)
        ! ...

        function read_from_file(this, filename) result(res)
            class(dXdrz_db_type), intent(inout) :: this
            character(*), intent(in) :: filename
            integer :: res

            character(len=256) :: line
            integer :: file_unit, read_pos
            logical :: reading
            type(dMdrz_type) :: dMdrz
            type(dMdrz_type), dimension(:), allocatable :: dMdrz_temp
            type(dAdrz_type) :: dAdrz
            type(dAdrz_type), dimension(:), allocatable :: dAdrz_temp

            res = 0

            open(newunit = file_unit, file = trim(adjustl(filename)), iostat = res)
            if(res .ne. 0) return
            
            reading = .false.
            read_pos = 0
            
            do while(res == 0)
                
                read(file_unit, "(A)", iostat = res) line
                if(res .ne. 0) cycle
                line = trim(adjustl(line))

                if(index(line, '#') == 1) cycle !found a comment
                if(index(line, ' ') == 1) cycle !found an empty line
                
                if(.not. reading) then
                    if(index(line, '$PARAMETERS') .ne. 1) then
                        res = -3
                        close(file_unit)
                        write(error_unit,"(A)") "DXDRZ_DB ERROR: invalid line: "//trim(line)
                        return !not a comment, not an empty line, not a token - error
                    end if
                    reading = .true.
                    cycle
                else
                    if(index(line, '$') == 1) then
                        res = -4
                        close(file_unit)
                        write(error_unit,"(A)") "DXDRZ_DB ERROR: unexpected token line: "//trim(line)
                        return !token in the middle of reading
                    end if
                    
                    if(read_pos == 0) then !read dMdrz
                        read(line,*, iostat = res) dMdrz%elements(:)
                        if(res .ne. 0) cycle
                        
                    else !read dAdrz
                        read(line,*, iostat = res) dAdrz%elements(read_pos, :)
                        if(res .ne. 0) cycle
                    end if

                    read_pos = read_pos + 1
                    
                    if(read_pos .gt. 3) then

                        !reallocate, append
                        if(allocated(dMdrz_temp)) deallocate(dMdrz_temp)
                        allocate(dMdrz_temp(this%count+1))
                        dMdrz_temp(1:this%count) = this%dMdrz_record(:)
                        dMdrz_temp(this%count+1) = dMdrz
                        call move_alloc(from = dMdrz_temp, to = this%dMdrz_record)

                        if(allocated(dAdrz_temp)) deallocate(dAdrz_temp)
                        allocate(dAdrz_temp(this%count+1))
                        dAdrz_temp(1:this%count) = this%dAdrz_record(:)
                        dAdrz_temp(this%count+1) = dAdrz
                        call move_alloc(from = dAdrz_temp, to = this%dAdrz_record)

                        this%count = this%count + 1
                        read_pos = 0
                        reading = .false.
                    end if
                    
                end if
            end do
            
            close(file_unit)

            if(reading) then
                res = -3
                write(error_unit,"(A)") "DXDRZ_DB ERROR: incomplete set of parameters"
                return !did not read the whole parameters list
            end if

            if(res == IOSTAT_END) then
                res = 0 !everything OK
                return
            end if

            write(error_unit,"(A)") "DXDRZ_DB ERROR: unexpected error"
        
        end function read_from_file
    
end module