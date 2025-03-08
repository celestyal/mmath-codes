program flux2d_gen_fields
    implicit none

    !
    ! File saving parameters
    !

    ! Only disable if you have already generated a profile and are looking to overwrite part of the profile
    logical, parameter :: save_ax               = .true.
    logical, parameter :: save_ay               = .true.
    logical, parameter :: save_vx               = .true.
    logical, parameter :: save_vy               = .true.
    logical, parameter :: save_grid             = .true.

    ! Directory to write output files to
    character(len=128), parameter :: directory  = "./exports/example/" ! Must have / at end
    
    !
    ! Grid parameters
    !

    ! Number of cells (higher -> more computational resources used)
    integer, parameter :: nx                    = 500
    integer, parameter :: ny                    = 500

    ! Dimensionless x and y minimum and maximum values 
    ! y=0 represents the equator
    ! x=0 has no physical importance
    ! Unit length represents the length scale, that is passed separately to the solver.
    real, parameter :: xmin                     = -0.5
    real, parameter :: xmax                     =  0.5
    real, parameter :: ymin                     = -0.5
    real, parameter :: ymax                     =  0.5

    !
    ! Magnetic vector potential configuration
    !

    ! 0 -> custom configuration. Modify custom_a to create the configuration.
    ! 1 -> current sheet
    ! 2 -> constant b
    ! 3 -> single bipole
    integer, parameter :: magnetic_profile      = 1

    ! Plasma velocity profile
    logical, parameter :: custom_flow           = .false. !Manually edit "custom_flow()"
    logical, parameter :: differential_rotation = .false.
    real, parameter :: diffrot_speed            = 0.0 !Speed at the equator of the star
    
    logical, parameter :: meridional_flow       = .false.
    real, parameter :: meridional_flow_speed    = 0.0
    
    logical, parameter :: convergent_flow       = .false.
    real, parameter :: convergent_x_coord       = 0.0
    real, parameter :: convergence_speed        = 0.0
    
    !!! End of modifiable parameters

    !!! Do not edit
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
            vx(:,i) = vx(:,i) + convergence_speed
        end do
        do concurrent (i = 1:nx+1, xcorn(i) .gt. convergent_x_coord)
            vx(:,i) = vx(:,i) - convergence_speed
        end do
    end subroutine convflow


    ! Differential rotation
    subroutine diffrot()
        real, parameter :: omega=2.0
        integer :: i

        do concurrent (i=1:nx+1)
            vx(:,i) = -omega*xrib(i)
        end do
    end subroutine diffrot


    ! Meridional (poleward) flow
    subroutine merflow()
        ! To be implemented
    end subroutine merflow


    !
    ! Do not edit
    ! Routines used to provide coordinates that can be used to generate
    ! the velocity and magnetic vector potential
    !
    subroutine coordinates()
        integer :: i
        real :: dx, dy

        dx = (xmax-xmin)/nx
        dy = (ymax-ymin)/ny

        do concurrent (i=1:nx+1)
            xcorn(i) = xmin + (dx * (i-1))
        end do
        do concurrent (i=1:ny+1)
            ycorn(i) = ymin + (dy * (i-1))
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
        real, parameter :: rho=0.20, b0=10.0
        integer :: i, j
        real, allocatable :: epsilon(:,:)

        allocate(epsilon(nx+1, ny))

        ax = 0.0
        do concurrent (i=1:nx+1, j=1:ny)
            epsilon(i,j) = (0.5 * (xcorn(i)**2) + yrib(j)**2) / (rho**2)
        end do

        ay = b0*rho*exp(0.5-epsilon)
    end subroutine single_bipole


    !
    ! Do not edit
    ! This code generates everything and saves the needed data used by the
    ! solver
    !
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
        !Custom
        case (0)
            call your_a()
        !Current Sheet
        case (1)
            call csheet()
        !Ax=-0.5*B0*y, Ay=0.5*B0*x (gives constant B)
        case (2)
            call const_b()
        !Single bipole
        case (3)
            call single_bipole()
        end select

        !Create output directory
        call execute_command_line('mkdir -p '//directory)

        !Write to file
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
end program flux2d_gen_fields
