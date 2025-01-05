program flux2d
    use flux2d_grid
    use flux2d_io
    use flux2d_ic
    use flux2d_solv
    use flux2d_types
    use flux2d_vel
    implicit none

    type(debug) :: d
    type(grid) :: g
    type(profile) :: p
    type(plane_component) :: a, v
    type(coordinates) :: ax_coord, ay_coord, grid_coord, b_coord
    real :: dx, dy, t1
    integer :: i
    character(512) :: id

    ! read the id in
    id = read_id()

    ! read the namelists in
    call read_profile(p, id)
    call read_debug(d, id)
    call read_grid(g, id)

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
    else
        ! need to read in the existing profile for a
    end if

    ! number of saves will be changed to the number of checkpoints. An additional start and end
    ! save will occur in addition to the checkpoints. Saving is currently not implemented.

    ! checks will need to made to ensure that there aren't more checkpoints than number of steps
    ! in the simulation.

    ! calculate save checkpoints
    t1 = (1/real(d%number_of_saves+1))*p%simulation_end_time 
    do i=1, d%number_of_saves
        ! solve
        call solv_ind(a%x, a%y, v%x, v%y, p%diffusivity, dx, dy, t1, p%boundary_condition, d%numerical_solver)
        print *, "checkpoint ", i, " of ", d%number_of_saves, " completed (no save was actually made as it its implementation is being rewritten)"
    end do
end program flux2d
