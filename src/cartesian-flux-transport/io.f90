!+
! :Description:
!   This module implements the parsing of command line arguments and the file
!   read and write routines. The main purpose of the module is to
!   - Read in the initial profile
!   - Enable the user to select simulation parameters when the program is
!     called
!   - Write saves to file
!   The program uses stream I/O to efficiently read and write unformatted data
!   to file. It is assumed that the program is compiled with a default integer
!   and real precision of 8 bytes (double precision). The program will not
!   run if it is compiled with different precisions.
!
!   The simulation file is structured as follows:
!     Preamble: version, header length
!     Header (each variable in the namelist is written separately)
!     vx
!     vy
!     t=t1 (initial save at t=0)
!     ax1
!     ay1
!     t=t2
!     ax2
!     ay2
!     ...
!     aysaves
!-

module cft_io
    implicit none

    !+ Header written to file
    type :: header ! 104 bytes
        integer :: solv, bc, nx, ny, sav
        real :: d, t, lx, ly, xmin, xmax, ymin, ymax
    end type header

    !+ Arguments used across the program are an extension of the header
    !  variables
    !-
    type, extends(header) :: args
        character(len=:), allocatable :: dir, file
        real :: freq
    end type args
contains
    subroutine reprof(ax, ay, vx, vy, t, arg)
        integer :: u
        real, allocatable, dimension(:,:) :: ax, ay, vx, vy
        logical :: exists
        type(args), intent(inout) :: arg
        real, intent(out) :: t

        if (arg%file .eq. '') then
            !+ Generate new simulation file if no existing file was supplied
            call newsim(ax, ay, vx, vy, arg)
        else
            !+ Load existing simulation
            call rehead(arg)
            call resav(arg%file, ax, ay, vx, vy, t=t, sav=-1)
        end if
    end subroutine reprof


    subroutine savsim(file, ax, ay, t)
        character(len=*), intent(in) :: file
        real, intent(in), dimension(:,:) :: ax, ay
        real, intent(in) :: t
        real :: told
        type(args) :: arg
        integer :: u, s

        !+ Chcek file existence
        call checkex(file)
        inquire(file=file, size=s)
        
        arg%file = file
        call rehead(arg)
        if (t .ne. 0) then
            !+ If not initial save, read previous save
            call resav(arg%file, t=told, sav=-1)
        else
            told = 0
        end if

        !+ Append save data to the simulation save file
        open(newunit=u, file=file, access='stream', status='old')
            write(unit=u, pos=s+1) t+told, ax, ay
        close(unit=u)

        !+ Update header with incremented save count
        arg%sav = arg%sav + 1
        call wrhead(file, arg)
    end subroutine savsim


    subroutine wrhead(file, arg)
        character(len=*), intent(in) :: file
        type(args) :: arg
        integer :: u
        
        !+ Check file exists
        call checkex(file)

        open(newunit=u, file=file, access='stream', status='old')
            !+ Write new/updated header to file
            write(u, pos=17) arg%solv, arg%bc, arg%nx, arg%ny, arg%sav, arg%d, arg%t, arg%lx, arg%ly, arg%xmin, arg%xmax, arg%ymin, arg%ymax
        close(unit=u)
    end subroutine wrhead


    subroutine wrpre(file)
        character(len=*), intent(in) :: file
        integer :: u, p
        integer, parameter :: version=1, hlen=104 

        !+ Write new/updated preamble to file
        open(newunit=u, file=file, access='stream', status='new')
            write(u, pos=1) version, hlen
        close(unit=u)
    end subroutine wrpre


    subroutine newsim(ax, ay, vx, vy, arg)
        real, allocatable, dimension(:,:), intent(out) :: ax, ay, vx, vy
        type(args), intent(inout) :: arg
        integer :: u
        
        !+ Load information about the grid from file
        call readgr(trim(arg%dir)//'grid_info.dat', arg%xmin, arg%xmax, arg%ymin, arg%ymax, arg%nx, arg%ny)

        !+ Allocate scalar fields
        allocate(ax(arg%nx, arg%ny+1))
        allocate(ay(arg%nx+1, arg%ny))
        allocate(vx(arg%nx+1, arg%ny+1))
        allocate(vy(arg%nx+1, arg%ny+1))

        !+ Read in scalar fields
        call rearr(trim(arg%dir)//'ax.dat', ax)
        call rearr(trim(arg%dir)//'ay.dat', ay)
        call rearr(trim(arg%dir)//'vx.dat', vx)
        call rearr(trim(arg%dir)//'vy.dat', vy)

        arg%file = trim(arg%dir)//'simulation.bin'

        !+ Write file preamble
        call wrpre(arg%file)

        !+ Write file header
        call wrhead(arg%file, arg)
        
        !+ Write velocities after header
        call wrvel(arg%file, vx, vy)

        !+ Make initial save (t=0)
        call savsim(arg%file, ax, ay, 0.0)
    end subroutine newsim


    subroutine wrvel(file, vx, vy)
        character(len=*), intent(in) :: file
        real, intent(in), dimension(:,:) :: vx, vy
        integer, parameter :: hlen=104
        integer :: u

        !+ Check file exists
        call checkex(file)

        !+ More robust checks needed here

        open(newunit=u, file=file, access='stream', status='old')
            !+ Write velocities to file after header
            write(u, pos=17+hlen) vx, vy
        close(unit=u)
    end subroutine wrvel


    subroutine rehead(arg)
        type(args), intent(inout) :: arg
        integer, parameter :: bypre=16
        integer :: u, v, l
        
        !+ Check file exists
        call checkex(arg%file)

        open(newunit=u, file=arg%file, status='old', action='read', access='stream', form='unformatted')
            !+
            ! Read preamble and then the header.
            ! The preamble is not currently used to deduce how to read the header ! but should be implemented later on.
            !-
            read(u, pos=1) v, l
            read(u, pos=bypre+1) arg%solv, arg%bc, arg%nx, arg%ny, arg%sav, arg%d, arg%t, arg%lx, arg%ly, arg%xmin, arg%xmax, arg%ymin, arg%ymax
        close(unit=u)
    end subroutine rehead


    subroutine revel(file, vx, vy)
        character(len=*), intent(in) :: file
        real, allocatable :: vx(:,:), vy(:,:)
        integer :: version, hlen
        integer :: u, offset

        open(newunit=u, file=file, status='old', action='read', access='stream', form='unformatted')
            read(u, pos=1) version, hlen
        
            offset = 16 + hlen + 1
            read(u, pos=offset) vx, vy
        close(u)
    end subroutine revel


    subroutine resav(file, ax, ay, vx, vy, t, sav)
        integer :: u, i, byax, byay, byv, byhead, byt
        real, allocatable, dimension(:,:), intent(out), optional :: ax, ay, vx, vy
        integer :: byf
        type(args) :: arg
        integer, parameter :: bypre=16, hlen=104 ! size of preamble
        integer, optional :: sav
        integer :: sav2
        integer :: offset
        character(len=*), intent(in) :: file
        real, optional :: t
        integer :: temp

        arg%file = file
        call rehead(arg)

        !+ Determine offsets needed to read file properly
        byax = arg%nx*(arg%ny+1)*8
        byay = (arg%nx+1)*(arg%ny)*8
        byv = (arg%nx+1)*(arg%ny+1)*8
        byhead = hlen
        byt = 8

        open(newunit=u, file=arg%file, status='old', action='read', access='stream', form='unformatted')
            !+ Read velocities
            offset = bypre + byhead + 1
            if (present(vx)) then
                allocate(vx(arg%nx+1, arg%ny+1))
                read(u, pos=offset) vx
            end if

            offset = offset + byv 
            if (present(vy)) then
                allocate(vy(arg%nx+1, arg%ny+1))
                read(u, pos=offset) vy
            end if

            !+ Read time, ax, ay of the given save
            offset = offset + byv
            if (present(sav)) then
                sav2 = sav
                if (sav .le. arg%sav) then
                    !+
                    ! If sav < 0, then access back of list.
                    ! -1 is the final element in the list.
                    ! -2 is the penultimate element, and so on.
                    !-
                    if (sav .lt. 0) then
                        sav2 = arg%sav + sav + 1
                        if (sav2 .gt. arg%sav) then
                            stop 'Attempting to read a save not present in the file, quitting.'
                        end if
                    end if
                    offset = offset + (sav2-1)*(byt+byax+byay)
                    !+ Read time
                    if (present(t)) then
                        read(u, pos=offset) t
                    end if
                    !+ Read ax
                    inquire(file=file, size=temp)
                    offset = offset + byt
                    if (present(ax)) then
                        allocate(ax(arg%nx, arg%ny+1))
                        read(u, pos=offset) ax
                    end if
                    !+ Read ay
                    offset = offset + byax
                    if (present(ay)) then
                        allocate(ay(arg%nx+1, arg%ny))
                        read(u, pos=offset) ay
                    end if
                end if
            end if
        close(unit=u)
    end subroutine resav


    subroutine rearr(file, data)
        character(len=*), intent(in) :: file
        real, allocatable, intent(inout) :: data(:,:)
        integer :: u

        !+ Check file exists
        call checkex(file)

        open(newunit=u, file=file, status='old', action='read', form='unformatted')
            read(u) data
        close(unit=u)
    end subroutine rearr


    subroutine readgr(file, xmin, xmax, ymin, ymax, nx, ny)
        character(len=*), intent(in) :: file
        real, intent(out) :: xmin, xmax, ymin, ymax
        integer, intent(out) :: nx, ny
        integer :: u

        open(newunit=u, file=file, status='old', action='read', form='unformatted')
            read(u) xmin, xmax, ymin, ymax, nx, ny
        close(unit=u)
    end subroutine readgr


    subroutine checkex(file)
        logical :: ex
        character(len=*), intent(in) :: file
        
        inquire(file=file, exist=ex)
        if (ex .eqv. .false.) then
            print *, file, ' does not exist. Cannot proceed.'
            stop
        end if
    end subroutine checkex


    function parse_args() result(arg)
        type(args) :: arg
        character(len=256), allocatable, dimension(:) :: arglist
        integer :: i, solver, nargs
        character(len=512) :: b

        !+ Check if arguments were supplied and print help dialogue if not
        nargs = command_argument_count()
        if (nargs .eq. 0) then
            stop
        endif

        allocate(arglist(nargs))

        do i=1, nargs
            call get_command_argument(i, arglist(i))
        end do

        !+
        ! Default arguments
        !-
        arg%solv = 2 ! Heun solver by default
        arg%freq = 0 ! Disable saving by default
        arg%file = ''
        arg%dir = ''
        arg%sav = 0
        arg%t = -1

        !+
        ! These will throw an error if not changed.
        ! This is intentional, to prevent simulations launching
        ! Without critical parameters being set
        !-
        arg%bc = -1
        arg%d = -1
        arg%t = -1
        arg%lx = -1
        arg%ly = -1

        !+ Parse command line arguments
        i=1
        do while (i .le. nargs)
            select case (arglist(i))
            case ('-s', '--freq')
                read(arglist(i+1), *) arg%freq
                i = i+1
            case ('--euler')
                arg%solv = 1
            case ('--heun')
                arg%solv = 2
            case ('--rk4')
                arg%solv = 3
            case ('-e', '--existing')
                arg%file = arglist(i+1)
                i = i+1
            case ('-n', '--new')
                arg%dir = arglist(i+1)
                i = i+1
            case ('-x', '--lx')
                read(arglist(i+1), *) arg%lx
                i = i+1
            case ('-y', '--ly')
                read(arglist(i+1), *) arg%ly
                i = i+1
            case ('-D', '--diffusivity', '-d')
                read(arglist(i+1), *) arg%d
                i = i+1
            case ('-b', '--boundary-condition', '--bc')
                b = arglist(i+1)
                select case (b)
                case ('periodic')
                    arg%bc = 1
                case ('infinite')
                    arg%bc = 2
                case ('open')
                    arg%bc = 3
                case default
                    print *, 'Not a valid boundary condition. Options are periodic, infinite, or open.'
                    stop
                end select
                i = i+1
            case ('-t', '--time')
                read(arglist(i+1), *) arg%t
                i = i + 1
            case default
                print '(a,a,/)', 'Ignoring unrecognized argument: ', arglist(i)
            end select
            i=i+1
        end do

        if (arg%file .eq. '') then
            if (arg%dir .eq. '') then
                print *, 'No existing or new simulation data provided. Provide -e "<PATH TO simulation binary FILE>" to load an existing simulation, or -n "<PATH TO SETUP DIRECTORY>" to start a new simulation.'
                stop
            end if
            if (arg%bc .eq. -1) then
                print *, "No boundary condition provided"
                stop
            end if
            
            if (arg%freq .eq. -1) then
                print *, "No save frequency (units of time) provided. Use -s 0 to disable checkpoints or -s <time> to select a frequency."
                stop
            end if

            if (arg%t .eq. -1) then
                print *, "No simulation time provided. Use -t <time> to set a time"
                stop
            end if
            
            if (arg%d .eq. -1) then
                print *, "No value for diffusivity provided. Use -d 0 to disable diffusion."
                stop
            end if
            
            if (arg%lx .eq. -1) then
                print *, "No x length scale provided. Use -x <length scale>."
                stop
            end if

            if (arg%ly .eq. -1) then
                print *, "No y length scale provided. Use -y <length scale>."
                stop
            end if
        end if
    end function parse_args
end module cft_io
