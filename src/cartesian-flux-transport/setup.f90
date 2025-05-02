program cft_setup
    implicit none
    !+
    ! Save the generated profiles to files
    ! Each of these are needed to set up the simulation
    !-
    logical, parameter :: save_ax               = .true.
    logical, parameter :: save_ay               = .true.
    logical, parameter :: save_vx               = .true.
    logical, parameter :: save_vy               = .true.
    logical, parameter :: save_grid             = .true.

    !+ Directory to write output files to, "./" is the working directory
    character(len=128), parameter :: directory  = "./"

    !+ Number of grid cells
    integer, parameter :: nx                    = 720
    integer, parameter :: ny                    = 360

    !+
    ! Minimum and maximum x and y values in the non-dimensionalised domain
    ! The length scale which is supplied to the solver is one unit of length
    ! when non-dimensionalised.
    !-
    real, parameter :: xmin                     = -2.0
    real, parameter :: xmax                     =  2.0
    real, parameter :: ymin                     = -1.0
    real, parameter :: ymax                     =  1.0

    !+
    ! Magnetic field initial condition
    ! Set:
    !  0 for a Custom configuration.
    !  (Modify your_a() to generate the configuration)
    !
    !  1 for a current sheet
    !  2 for constant b
    !  3 for a single bipole
    !-
    integer, parameter :: magnetic_profile      = 3

    !+ A custom flow can be superposed with the other flows
    logical, parameter :: custom_flow           = .false.

    !+
    ! Differential rotation is assumed to be taken in a synodic frame
    ! The reference point represents the latitude on the Sun where the rate of
    ! rotation mathces the rate that the frame of reference moves
    !-
    logical, parameter :: differential_rotation = .false.
    real, parameter :: diffrot_coefficient      = 0.2
    real, parameter :: reference_point          = 0.0

    ! Meridional flow is taken to be contant across the domain.
    logical, parameter :: meridional_flow       = .false.
    real, parameter :: meridional_flow_speed    = 0.011
    
    ! Convergent flow transports magnetic field towards a particular coordinate.
    logical, parameter :: convergent_flow       = .false.
    real, parameter :: convergent_x_coord       = 0.0
    real, parameter :: convergence_speed        = 0.0
    
    !+
    ! No need to edit beyond this line unless implementing a custom magnetic
    ! field or flow.
    !-

    real, dimension(nx, ny+1) :: ax
    real, dimension(nx) :: xrib
    real, dimension(ny) :: yrib
    real, dimension(nx+1) :: xcorn
    real, dimension(ny+1) :: ycorn
    real, dimension(nx+1, ny) :: ay
    real, dimension(nx+1, ny+1) :: vx, vy
    integer :: u
    
    ! Initialise quantities to 0
    ax = 0.0
    ay = 0.0
    vx = 0.0
    vy = 0.0

    ! Call code to generate and save quantities
    call generator()
    !!! End of do not edit
contains
    !!! Modify your_a and your_flow if using custom routines
    subroutine your_a()
        ! Insert your code to populate for ax, ay
    end subroutine your_a


    subroutine your_flow()
        ! Insert your code to populate vx, vy
    end subroutine your_flow
    !!! End of custom profiles

    !
    ! Subroutines used to implement the effects of the different solar
    ! flows. These can be edited to change the way these flows are implemented
    !

    ! Convergent flow
    subroutine convflow()
        integer :: i

        do concurrent (i = 1:nx+1, xcorn(i) .lt. convergent_x_coord)
            vx(i,:) = vx(i,:) + convergence_speed
        end do
        do concurrent (i = 1:nx+1, xcorn(i) .gt. convergent_x_coord)
            vx(i,:) = vx(i,:) - convergence_speed
        end do
    end subroutine convflow


    ! Differential rotation
    subroutine diffrot()
        integer :: i

        !
        ! Taken to be negative at high latitudes, accounting for the Earth's orbit
        ! in a sidereal frame of reference.
        ! The stagnation point is the point where the reversal occurs (at mid latitudes)
        !
        do concurrent (i=1:ny+1)
            vx(:, i) = vx(:, i) - diffrot_coefficient*(ycorn(i)-reference_point)
        end do
    end subroutine diffrot


    ! Meridional (poleward) flow
    subroutine merflow()
        vy = vy + meridional_flow_speed
    end subroutine merflow


    ! Coordinates needed to help populate vector potential and velocities
    subroutine coordinates()
        integer :: i
        real :: dx, dy

        dx = (xmax-xmin)/nx
        dy = (ymax-ymin)/ny

        do concurrent (i=1:nx+1)
            xcorn(i) = xmin + dx*(i-1)
        end do
        do concurrent (i=1:ny+1)
            ycorn(i) = ymin + dy*(i-1)
        end do

        xrib = xcorn(2:) - dx*0.5
        yrib = ycorn(2:) - dy*0.5
    end subroutine coordinates


    ! Continuous current sheet
    subroutine csheet()
        integer :: i
        real, parameter :: b0 = 250.0
        
        ax = 0.0
        do concurrent (i=1:nx+1)
            ay(i,:) = (10/b0)*log(cosh(b0 * xcorn(i)))
        end do
    end subroutine csheet


    ! Consant B (for testing)
    subroutine const_b()
    integer :: i
        real, parameter :: b0 = 10.0

        do concurrent (i=1:nx+1)
            ay(i,:) = 0.5*b0*xcorn(i)
        end do
        do concurrent (i=1:ny+1)
            ax(:,i) = -0.5*b0*ycorn(i)
        end do
    end subroutine const_b


    ! Single bipole
    subroutine single_bipole()
        real, parameter :: rho=0.40, b0=-50.0, tilt=(4.0*ATAN(1.0)/180.0)*25.0, x0=0, y0=0
        integer :: i, j
        real :: epsilon(nx+1, ny), xprime(nx+1, ny), yprime(nx+1, ny)

        ax = 0.0

        ! Apply transform (to produce tilt)
        do concurrent (i=1:nx+1, j=1:ny)
            xprime(i, j) = (xcorn(i)-x0)*cos(tilt) - (yrib(j)-y0)*sin(tilt)
            yprime(i, j) = (xcorn(i)-x0)*sin(tilt) + (yrib(j)-y0)*cos(tilt)
        end do

        epsilon = (0.5 * (xprime**2) + yprime**2) / (rho**2)
        ay = b0*rho*exp(0.5-epsilon)
    end subroutine single_bipole


    ! The generator subroutine applies the configuration parameters
    subroutine generator()
        call coordinates()

        if (custom_flow .eqv. .true.) then
            call your_flow()
        end if
        if (differential_rotation .eqv. .true.) then
            call diffrot()
        end if
        if (meridional_flow .eqv. .true.) then
            call merflow()
        end if
        if (convergent_flow .eqv. .true.) then
            call convflow()
        end if

        select case (magnetic_profile)
        ! Custom
        case (0)
            call your_a()
        ! Current Sheet
        case (1)
            call csheet()
        ! Ax=-0.5*B0*y, Ay=0.5*B0*x (gives constant B)
        case (2)
            call const_b()
        ! Single bipole
        case (3)
            call single_bipole()
        end select

        ! Create output directory if not exists (only tested with GNU coreutils)
        call execute_command_line('mkdir -p '//directory)

        ! Write to file
        if (save_ax .eqv. .true.) then
            open(newunit=u, file=trim(directory)//"ax.dat", form='unformatted', status='new', action='write')
                write(u) ax
            close(u)
        endif

        if (save_ay .eqv. .true.) then
            open(newunit=u, file=trim(directory)//"ay.dat", form='unformatted', status='new', action='write')
                write(u) ay
            close(u)
        endif

        if (save_vx .eqv. .true.) then
            open(newunit=u, file=trim(directory)//"vx.dat", form='unformatted', status='new', action='write')
                write(u) vx
            close(u)
        endif

        if (save_vy .eqv. .true.) then
            open(newunit=u, file=trim(directory)//"vy.dat", form='unformatted', status='new', action='write')
                write(u) vy
            close(u)
        endif

        if (save_grid .eqv. .true.) then
            open(newunit=u, file=trim(directory)//"grid_info.dat", form='unformatted', status='new', action='write')
                write(u) xmin, xmax, ymin, ymax, nx, ny
            close(u)
        endif
    end subroutine generator
end program cft_setup
