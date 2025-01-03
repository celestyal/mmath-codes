program flux_2d
    use mag_grid
    use file_io
    use vector_potential
    use numerical_solver
    use types
    use velocity
    implicit none

    type(debug) :: d
    type(grid) :: g
    type(profile) :: p
    type(plane_component) :: a, v
    type(coordinates) :: ax_coord, ay_coord, grid_coord, b_coord
    real :: dt, dx, dy
    character(128) :: file, int_char
    character(512) :: identifier

    ! read the identifier in
    identifier =  read_identifier()

    ! read the namelists in
    call read_profile(p, identifier)
    call read_debug(d, identifier)
    call read_grid(g, identifier)

    ! derive dx and dy
    dx = (g%xmax-g%xmin)/g%nx
    dy = (g%ymax-g%ymin)/g%ny

    ! allocate stuff using grid variables
    allocate(a%x(g%nx, g%ny+1))
    allocate(a%y(g%nx+1, g%ny))
    allocate(ax_coord%x(g%nx))
    allocate(ax_coord%y(g%ny+1))
    allocate(ay_coord%x(g%nx+1))
    allocate(ay_coord%y(g%ny))
    allocate(grid_coord%x(g%nx+1))
    allocate(grid_coord%y(g%ny+1))
    allocate(b_coord%x(g%nx+2))
    allocate(b_coord%y(g%ny+2))

    call populate_grid_coordinates(ax_coord%x, ax_coord%y, ay_coord%x, &
        ay_coord%y, grid_coord%x, grid_coord%y, b_coord%x, b_coord%y, &
        g%xmin, g%ymin, dx, dy, g%nx, g%ny)

    allocate(v%x(g%nx+1, g%ny+1))
    allocate(v%y(g%nx+1, g%ny+1))
    v%x = 0.0
    v%y = 0.0

    ! populate velocities
    call apply_convergent_flow(v%x, p%convergence_speed, p%convergence_point, grid_coord%x, g%nx)
    call apply_differential_flow(v%x, p%differential_speed, grid_coord%y, g%nx)
    call apply_meridional_flow(v%y, p%meridional_speed)

    ! scale velocities to be units per second
    call apply_velocity_scaling(p%diffusivity, v%x, v%y, p%length_scale)

    ! find timestep that satisfies cfl condition
    call find_dt(dt, v%x, v%y, dx, dy, p%diffusivity)
    p%steps = nint( (p%simulation_end_time-p%simulation_start_time)/dt )

    ! initialise quantity arrays
    call init_farrays(g%nx, g%ny)

    if (d%new_simulation) then
        d%new_simulation = .false.
        if (p%potential_configuration .eq. 0) then
            !we will load a custom configuration
        else if (p%potential_configuration .eq. -1) then
            !we will read a configuration for bz and then compute ax and ay from that
        else if (p%potential_configuration .eq. 1) then !load current sheet configuration
            call current_sheet(a%x, a%y, ay_coord%x, g%nx)
        else if (p%potential_configuration .eq. 2) then !load advection test
            call advection_test(a%x, a%y, ay_coord%x, ax_coord%y, g%nx, g%ny)
        else if (p%potential_configuration .eq. 3) then !load bipole with no tilt
            call bipole_no_tilt(a%x, a%y, ay_coord%x, ay_coord%y, g%nx, g%ny)
        end if

        ! export initial snapshot of the system before evolution
        if (d%number_of_saves .gt. 0) then
            call save_snapshot(.true., identifier, a%x, a%y, dx, dy, g%nx, g%ny, vx=v%x, vy=v%y)
            call export_for_idl(identifier, g)
            call save_timestamp(identifier, p%simulation_start_time)
        end if
    else
        ! need to read in the existing profile for a
    end if

    ! use numerical solver
    select case (d%numerical_solver)
    ! euler
    case (1)
        call euler(a%x, a%y, v%x, v%y, p%diffusivity, p%boundary_condition, dx, dy, dt, p%steps, g%nx, g%ny, p%simulation_start_time, d%number_of_saves, identifier)

    ! heun
    case (2)
        call heun(a%x, a%y, v%x, v%y, p%diffusivity, p%boundary_condition, dx, dy, dt, p%steps, g%nx, g%ny, p%simulation_start_time, d%number_of_saves, identifier)

    ! runge-kutta 4
    case (3)
        call rk4(a%x, a%y, v%x, v%y, p%diffusivity, p%boundary_condition, dx, dy, dt, p%steps, g%nx, g%ny, p%simulation_start_time, d%number_of_saves, identifier)

    case default
        error stop
    end select

contains

    subroutine find_dt(dt, vx, vy, dx, dy, diffusivity)
        real :: mvx, mvy
        real, dimension(3) :: t
        real, intent(in) :: vx(:,:), vy(:,:), dx, dy, diffusivity
        real, intent(out) :: dt
        integer :: i

        ! max velocities (km/s)
        mvx = maxval(abs(vx))
        mvy = maxval(abs(vy))

        ! times to cross each grid
        t(1) = dx / mvx
        t(2) = dy / mvy
        t(3) = (dx*dy) / abs(diffusivity)

        ! reference value to compare to
        dt = huge(1.)

        do concurrent (i = 1:3, (t(i) .gt. 0) .and. (t(i) .lt. dt))
        dt = t(i)
        end do

        dt = dt / 15.0
    end subroutine find_dt
end program flux_2d
