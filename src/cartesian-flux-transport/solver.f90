!+
! :Description:
!   This program solves Cartesian Flux Transport equations that were derived
!   as a Cartesian analogue to the standard spherical model. The program
!   expects the initial scalar fields (ax, ay, vx, vy) and grid info to be
!   provided through files. Other simulation parameters are passed to the
!   program using a command line input
!
!   The program:
!   - Parses the command line arguments and reads/creates a simulation,
!   - Selects a timestep that satisfies the CFL condition.
!   - Partitions the total simulation time into intervals to enable the
!     saving of the system at regular points during solving.
!-
program cft_solver
    use cft_io
    use cft_kernel
    implicit none

    integer :: i, iter
    character(128) :: id
    type(args) :: arg
    real, allocatable, dimension(:,:) :: ax, ay, vx, vy, ax0, ay0
    real :: tcurr, dt, multiplier, dx, dy, tnew

    !+ Read arguments supplied to the executable
    arg = parse_args()

    !+ Modify simulation time in file if -t is provided with an existing simulation
    if ((arg%t .ne. -1) .and. (arg%file .ne. '')) then
        tnew = arg%t
        call reprof(ax, ay, vx, vy, tcurr, arg)
        arg%t = tnew
    else
        !+ Load existing simulation or set up a new simulation
        call reprof(ax, ay, vx, vy, tcurr, arg)
    end if

    if (tcurr .lt. arg%t) then
        !+ Find a suitably small dt (necessary for convergence)
        dx = arg%lx*(arg%xmax-arg%xmin)/arg%nx
        dy = arg%ly*(arg%ymax-arg%ymin)/arg%ny
        dt = find_dt(vx, vy, dx, dy, arg%d)/5.0

        if (arg%freq .gt. 0) then
            !+ Select the lower of dt or the save frequency
            dt = min(dt, arg%freq)
            !+
            ! Adjust dt so that the save frequency is an exact multiple of dt.
            ! This is to force a save to occur only at multiples of freq.
            ! The final save will occur at the next save frequency after the
            ! given save time instead of the closest time using the original
            ! dt.
            !-
            multiplier = ceiling(arg%freq/dt)
            dt = arg%freq/multiplier
        else
            !+ Set freq to be the total remaining simulation time
            arg%freq = arg%t - tcurr
        end if

        !+ Read in the initial condition if a fixed boundary is used
        if (arg%bc .eq. 2) then
           allocate(ax0(arg%nx, arg%ny+1))
           allocate(ay0(arg%nx+1, arg%ny))
           call resav(arg%file, ax0, ay0, sav=1)
        end if
        
        !+
        ! Solve until reaching the next save point
        ! Further work: dynamic simulation parameters could be implemented
        ! by updating the profile during this loop.
        !-
        if (arg%bc .eq. 2) then
            do while (tcurr .lt. arg%t)
                ! Pass ax0 and ay0
                call solv_ind(ax, ay, vx, vy, arg%d, arg%lx, arg%ly, dt, arg%freq, arg%bc, arg%solv, arg%xmin, arg%xmax, arg%ymin, arg%ymax, ax0=ax0, ay0=ay0)
                call savsim(arg%file, ax, ay, arg%freq)
                tcurr = tcurr + arg%freq
                print *, "Save ", arg%sav, "made at t =", tcurr
            end do
        else
            do while (tcurr .lt. arg%t)
                ! Don't pass ax0 and ay0
                call solv_ind(ax, ay, vx, vy, arg%d, arg%lx, arg%ly, dt, arg%freq, arg%bc, arg%solv, arg%xmin, arg%xmax, arg%ymin, arg%ymax)
                call savsim(arg%file, ax, ay, arg%freq)
                tcurr = tcurr + arg%freq
                print *, "Save made at t =", tcurr
            end do
        end if
    else
        print *, "Simulation already completed."
    end if
end program cft_solver
