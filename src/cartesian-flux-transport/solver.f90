program cft_solver
    use cft_io
    use cft_kernel
    implicit none

    integer :: i, iter
    character(128) :: id
    type(args) :: arg
    real, allocatable, dimension(:,:) :: ax, ay, vx, vy
    real :: tcurr, dt, multiplier, dx, dy, tnew

    ! Read arguments supplied to the executable
    arg = parse_args()

    ! Modify simulation time in file if -t is provided with an existing simulation
    if ((arg%t .ne. -1) .and. (arg%file .ne. '')) then
        tnew = arg%t
        call reprof(ax, ay, vx, vy, tcurr, arg)
        arg%t = tnew
    else
        ! Load existing simulation or set up a new simulation
        call reprof(ax, ay, vx, vy, tcurr, arg)
    end if

    if (tcurr .lt. arg%t) then
        ! Find a suitably small dt (necessary for convergence)
        dx = arg%lx*(arg%xmax-arg%xmin)/arg%nx
        dy = arg%ly*(arg%ymax-arg%ymin)/arg%ny
        dt = find_dt(vx, vy, dx, dy, arg%d)

        if (arg%freq .gt. 0) then
            ! Select the lower of dt or the save frequency
            dt = min(dt, arg%freq)
            ! Adjust dt so that the save frequency is an exact multiple of dt.
            ! This is to force a save to occur only at multiples of freq.
            ! The final save will occur at the next save frequency after the given save time instead of the closest time using the
            ! original dt.
            multiplier = ceiling(arg%freq/dt)
            dt = arg%freq/multiplier
        else
            ! Set freq to be the total remaining simulation time
            arg%freq = arg%t - tcurr
        end if

        
        do while (tcurr .lt. arg%t)
            ! Solve until reaching the next save point
            call solv_ind(ax, ay, vx, vy, arg%d, arg%lx, arg%ly, dt, arg%freq, arg%bc, arg%solv, arg%xmin, arg%xmax, arg%ymin, arg%ymax)
            call savsim(arg%file, ax, ay, arg%freq)
            tcurr = tcurr + arg%freq
            print *, "Save ", arg%sav, "made at t =", tcurr
        end do
    else
        print *, "Simulation already completed."
    end if
end program cft_solver
