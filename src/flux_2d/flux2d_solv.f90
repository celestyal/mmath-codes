module flux2d_solv
    implicit none

    real, private, allocatable, dimension(:,:) :: b, jx, jy, b_corner, vxx, vxy, vxx_yrib, vxy_xrib
    integer, private :: nx, ny
contains
    subroutine solv_ind(ax, ay, vx, vy, d, dx, dy, t, bc, solver)
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy, t
        real, intent(inout) :: ax(:,:), ay(:,:)
        integer, intent(in) :: bc, solver
        real :: dt

        ! determine grid size using A to save passing unneeded variables in
        nx = size(ax, 1)
        ny = size(ay, 2)

        ! determine timestep
        dt = find_dt(vx, vy, dx, dy, d)

        ! initialise relevant arrays
        if ( allocated(b_corner) .eqv. .false. ) then
            allocate(b_corner(nx+1, ny+1))
            allocate(vxx(nx+1, ny+1))
            allocate(vxy(nx+1, ny+1))
            allocate(vxx_yrib(nx+1, ny))
            allocate(vxy_xrib(nx, ny+1))
            allocate(b(nx+2,ny+2))
            allocate(jx(nx,ny+1))
            allocate(jy(nx+1,ny))
        end if

        ! call desired numerical solver method
        select case (solver)
        case (1) 
            call euler(ax, ay, vx, vy, d, bc, dx, dy, dt, t)
        case (2)
            call heun(ax, ay, vx, vy, d, bc, dx, dy, dt, t)
        case (3)
            call rk4(ax, ay, vx, vy, d, bc, dx, dy, dt, t)
        case default
            call heun(ax, ay, vx, vy, d, bc, dx, dy, dt, t)
        end select
    end subroutine solv_ind


    subroutine f(ax, ay, vx, vy, d, bc, dx, dy, dax, day)
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy
        integer, intent(in) :: bc
        real, intent(in) :: ax(:,:), ay(:,:)
        real, intent(out) :: dax(:,:), day(:,:)

        ! Calculate b
        b(2:nx+1, 2:ny+1) = (ay(2:nx+1, 1:ny)-ay(1:nx, 1:ny)) / dx - (ax(1:nx, 2:ny+1)-ax(1:nx, 1:ny)) / dy
        
        ! Apply boundary condition to left and right boundaries
        select case (bc)
        case (1)
            ! periodic
            b(1,:) = b(nx+1,:)
            b(nx+2,:) = b(2,:)
        case (2)
            !infinite
            b(nx+2,:) = 10 !update this
            b(1,:) = -10
        case (3)
            b(nx+2,:) = b(nx+1,:)
            b(1,:) = b(2,:)
        end select

        ! apply a closed boundary condition at top and botton
        b(:,ny+2) = b(:,ny+1)
        b(:,1) = b(:,2)

        ! 4 point average for b at corner
        b_corner = (b(:nx+1, :ny+1) + b(:nx+1, 2:) + b(2:, :ny+1) + b(2:, 2:)) / 4.0

        ! on the corners
        vxx = b_corner * vx
        vxy = b_corner * vy

        ! on the locations of ax and ay
        vxx_yrib = (vxx(:, :ny) + vxx(:, 2:)) / 2.0
        vxy_xrib = (vxy(:nx, :) + vxy(2:, :)) / 2.0

        ! partial derivatives of jx and jy
        jx = (b(2:nx+1, 2:) - b(2:nx+1, :ny+1)) / dy
        jy =-(b(2:, 2:ny+1) - b(:nx+1, 2:ny+1)) / dx

        ! induction equation
        dax = vxy_xrib - d*jx
        day = -vxx_yrib - d*jy
    end subroutine f


    subroutine rk4(ax, ay, vx, vy, d, bc, dx, dy, dt, t)
        real, intent(inout) :: ax(:,:), ay(:,:)
        integer, intent(in) :: bc
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy, dt, t
        real, dimension(nx, ny+1) :: x1, x2, x3, x4
        real, dimension(nx+1, ny) :: y1, y2, y3, y4
        integer :: i

        do i=1, ceiling(t/dt)
            call f(ax, ay, vx, vy, d, bc, dx, dy, x1, y1)
            call f(ax + 0.5*dt*x1, ay + 0.5*dt*y1, vx, vy, d, bc, dx, dy, x2, y2)
            call f(ax + 0.5*dt*x2, ay + 0.5*dt*y2, vx, vy, d, bc, dx, dy, x3, y3)
            call f(ax + dt*x3, ay + dt*y3, vx, vy, d, bc, dx, dy, x4, y4)
            ax = ax + (dt/6)*(x1 + 2*x2 + 2*x3 + x4)
            ay = ay + (dt/6)*(y1 + 2*y2 + 2*y3 + y4)
        end do
    end subroutine rk4


    subroutine euler(ax, ay, vx, vy, d, bc, dx, dy, dt, t)
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy, dt, t
        real, intent(inout) :: ax(:,:), ay(:,:)
        integer, intent(in) :: bc
        real :: kx(nx, ny+1), ky(nx+1, ny)
        integer :: i

        do i=1, ceiling(t/dt)
            call f(ax, ay, vx, vy, d, bc, dx, dy, kx, ky)
            ax = ax + dt*kx
            ay = ay + dt*ky
        end do
        print *, t/dt
    end subroutine euler


    subroutine heun(ax, ay, vx, vy, d, bc, dx, dy, dt, t)
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy, dt, t
        real, intent(inout) :: ax(:,:), ay(:,:)
        integer, intent(in) :: bc
        real :: k1x(nx, ny+1), k1y(nx+1, ny), k2x(nx, ny+1), k2y(nx+1, ny)
        integer :: i

        do i=1, ceiling(t/dt)
            call f(ax, ay, vx, vy, d, bc, dx, dy, k1x, k1y)
            call f(ax + dt*k1x, ay + dt*k1y, vx, vy, d, bc, dx, dy, k2x, k2y)
            ax = ax + 0.5*dt*(k1x+k2x)
            ay = ay + 0.5*dt*(k1y+k2y)
        end do
    end subroutine heun


    pure function find_dt(vx, vy, dx, dy, d) result(dt)
        real :: mvx, mvy, dt
        real, dimension(3) :: t
        real, intent(in) :: vx(:,:), vy(:,:), dx, dy, d
        integer :: i

        ! max velocities (km/s)
        mvx = maxval(abs(vx))
        mvy = maxval(abs(vy))

        ! times to cross each grid
        t(1) = dx / mvx
        t(2) = dy / mvy
        t(3) = (dx*dy) / abs(d)

        ! reference value to compare to
        dt = huge(1.)

        do concurrent (i = 1:3, (t(i) .gt. 0) .and. (t(i) .lt. dt))
            dt = t(i)
        end do

        dt = dt / 10.0
    end function find_dt
end module flux2d_solv
