!+
! :Description:
!   This module implements the numerical methods which solve the Cartesian
!   Surface Flux Transport equations.
!-
module cft_kernel
    implicit none

    real, private, allocatable, dimension(:,:) :: b, jy, jx, b_corner, vxx, vxy, vxx_yrib, vxy_xrib, b0

    private :: f, rk4, euler, heun

    integer, private :: nx, ny
contains
    ! This is the subroutine that should be interfaced with to make use of this module
    subroutine solv_ind(ax, ay, vx_orig, vy_orig, d_orig, lx, ly, dt, t, bc, solver, xmin, xmax, ymin, ymax, ax0, ay0)
        real, intent(in) :: vx_orig(:,:), vy_orig(:,:), d_orig, lx, ly, dt, t, xmin, xmax, ymin, ymax
        real, intent(inout) :: ax(:,:), ay(:,:)
        integer, intent(in) :: bc, solver
        real, optional, intent(in) :: ax0(:,:), ay0(:,:)
        real :: d, dx, dy
        real, allocatable, dimension(:,:) :: vx, vy

        ! determine grid size using A to save passing unneeded variables in
        nx = size(ax, 1)
        ny = size(ay, 2)

        allocate(vx(nx+1, ny+1))
        allocate(vy(nx+1, ny+1))

        !+ Scale velocities to units/second
        vx = vx_orig/lx
        vy = vy_orig/ly

        dx = (xmax-xmin)/nx
        dy = (ymax-ymin)/ny

        ! Scale diffusion to units^2/second
        d = d_orig/(lx*ly)

        !+ Allocate array memory
        allocate(b_corner(nx+1, ny+1))
        allocate(vxx(nx+1, ny+1))
        allocate(vxy(nx+1, ny+1))
        allocate(vxx_yrib(nx+1, ny))
        allocate(vxy_xrib(nx, ny+1))
        allocate(b(nx+2, ny+2))
        allocate(jx(nx, ny+1))
        allocate(jy(nx+1, ny))
        if (present(ax0) .and. present(ay0)) then
            allocate(b0(nx+2, ny+2))
        end if

        !+ Call a numerical solver
        if (present(ax0) .and. present(ay0)) then
            select case (solver)
            case (1)
                call euler(ax, ay, vx, vy, d, bc, dx, dy, dt, t, ax0=ax0, ay0=ay0)
            case (2)
                call heun(ax, ay, vx, vy, d, bc, dx, dy, dt, t, ax0=ax0, ay0=ay0)
            case (3)
                call rk4(ax, ay, vx, vy, d, bc, dx, dy, dt, t, ax0=ax0, ay0=ay0)
            case default
                call heun(ax, ay, vx, vy, d, bc, dx, dy, dt, t, ax0=ax0, ay0=ay0)
            end select
        else
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
        end if
        
        !+
        ! Deallocate array memory (to prevent an attempt to reallocate them
        ! in subsequent calls
        !-
        deallocate(b_corner)
        deallocate(vxx)
        deallocate(vxy)
        deallocate(vxx_yrib)
        deallocate(vxy_xrib)
        deallocate(b)
        if (present(ax0) .and. present(ay0)) then
            deallocate(b0)
        end if
        deallocate(jx)
        deallocate(jy)
    end subroutine solv_ind

    !+ Calculates the right hand side of the model equations
    subroutine f(ax, ay, vx, vy, d, bc, dx, dy, dax, day, ax0, ay0)
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy
        integer, intent(in) :: bc
        real, optional, intent(in) :: ax0(:,:), ay0(:,:)
        real, intent(in) :: ax(:,:), ay(:,:)
        real, intent(out) :: dax(:,:), day(:,:)

        !+ Calculate B at cell centres
        b(2:nx+1, 2:ny+1) = (ay(2:nx+1, 1:ny)-ay(1:nx, 1:ny)) / dx - (ax(1:nx, 2:ny+1)-ax(1:nx, 1:ny)) / dy
        
        !+ Apply boundary condition
        select case (bc)
        !+ Periodic
        case (1)
            b(1,:) = b(nx+1,:)
            b(nx+2,:) = b(2,:)
        !+ Fixed
        case (2)
            if (present(ax0) .and. present(ay0)) then
                b0(2:nx+1, 2:ny+1) = (ay0(2:nx+1, 1:ny)-ay0(1:nx, 1:ny)) / dx - (ax0(1:nx, 2:ny+1)-ax0(1:nx, 1:ny)) / dy

                b(nx+2,:) = b0(nx+1,:)
                b(1,:) = b0(2,:)
            else
                print *, "No fixed boundary profile provided. Cannot proceed"
                stop
            endif

        !+ Open
        case (3)
            b(nx+2,:) = b(nx+1,:)
            b(1,:) = b(2,:)
        end select

        !+ Apply a closed boundary at top and botton
        b(:,ny+2) = b(:,ny+1)
        b(:,1) = b(:,2)

        !+ 4 point average to find magnetic field at cell corners
        b_corner = (b(:nx+1, :ny+1) + b(:nx+1, 2:) + b(2:, :ny+1) + b(2:, 2:)) / 4.0

        !+ Multiply velocity with magnetic field at the cell corners
        vxx = b_corner * vx
        vxy = b_corner * vy

        !+ 2 point average to find on the cell edge midpoints
        vxx_yrib = (vxx(:, :ny) + vxx(:, 2:)) / 2.0
        vxy_xrib = (vxy(:nx, :) + vxy(2:, :)) / 2.0

        !+ Partial derivatives of B (to find curl B)
        jx = (b(2:nx+1, 2:) - b(2:nx+1, :ny+1)) / dy
        jy =-(b(2:, 2:ny+1) - b(:nx+1, 2:ny+1)) / dx

        !+ RHS of induction equation
        dax = vxy_xrib - d*jx
        day = -vxx_yrib - d*jy
    end subroutine f


    !+ Chooses a value for dt that satisfied the CFL condition
    !- Note: this is necessary but not sufficient for convergence
    function find_dt(vx, vy, dx, dy, d) result(dt)
        real :: dt
        real, dimension(3) :: t
        real, intent(in) :: vx(:,:), vy(:,:), dx, dy, d
        integer :: i

        ! Times to to cross a cell
        t(1) = dx / maxval(abs(vx))
        t(2) = dy / maxval(abs(vy))
        t(3) = (dx*dy) / abs(d)

        ! Reference value to compare to
        dt = huge(1.)

        !+ Choose the smallest non-zero value
        do concurrent (i = 1:3, (t(i) .gt. 0) .and. (t(i) .lt. dt))
            dt = t(i)
        end do
    end function find_dt


    !+ Runge Kutta 4 solver
    subroutine rk4(ax, ay, vx, vy, d, bc, dx, dy, dt, t, ax0, ay0)
        real, intent(inout) :: ax(:,:), ay(:,:)
        integer, intent(in) :: bc
        real, optional, intent(in) :: ax0(:,:), ay0(:,:)
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy, dt, t
        real, dimension(nx, ny+1) :: x1, x2, x3, x4
        real, dimension(nx+1, ny) :: y1, y2, y3, y4
        integer :: i

        if (present(ax0) .and. present(ay0)) then
            !+ Don't try to pass b0
            do i=1, ceiling(t/dt)
                call f(ax, ay, vx, vy, d, bc, dx, dy, x1, y1, ax0=ax0, ay0=ay0)
                call f(ax + 0.5*dt*x1, ay + 0.5*dt*y1, vx, vy, d, bc, dx, dy, x2, y2, ax0=ax0, ay0=ay0)
                call f(ax + 0.5*dt*x2, ay + 0.5*dt*y2, vx, vy, d, bc, dx, dy, x3, y3, ax0=ax0, ay0=ay0)
                call f(ax + dt*x3, ay + dt*y3, vx, vy, d, bc, dx, dy, x4, y4, ax0=ax0, ay0=ay0)
                ax = ax + (dt/6)*(x1 + 2*x2 + 2*x3 + x4)
                ay = ay + (dt/6)*(y1 + 2*y2 + 2*y3 + y4)
            end do
        else
            do i=1, ceiling(t/dt)
                call f(ax, ay, vx, vy, d, bc, dx, dy, x1, y1)
                call f(ax + 0.5*dt*x1, ay + 0.5*dt*y1, vx, vy, d, bc, dx, dy, x2, y2)
                call f(ax + 0.5*dt*x2, ay + 0.5*dt*y2, vx, vy, d, bc, dx, dy, x3, y3)
                call f(ax + dt*x3, ay + dt*y3, vx, vy, d, bc, dx, dy, x4, y4)
                ax = ax + (dt/6)*(x1 + 2*x2 + 2*x3 + x4)
                ay = ay + (dt/6)*(y1 + 2*y2 + 2*y3 + y4)
            end do
        end if
    end subroutine rk4


    !+ Euler solver
    subroutine euler(ax, ay, vx, vy, d, bc, dx, dy, dt, t, ax0, ay0)
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy, dt, t
        real, intent(inout) :: ax(:,:), ay(:,:)
        integer, intent(in) :: bc
        real, optional, intent(in) :: ax0(:,:), ay0(:,:)
        real :: kx(nx, ny+1), ky(nx+1, ny)
        integer :: i

        if (present(ax0) .and. present(ay0)) then
            !+ Don't try to pass b0
            do i=1, ceiling(t/dt)
                call f(ax, ay, vx, vy, d, bc, dx, dy, kx, ky)
                ax = ax + dt*kx
                ay = ay + dt*ky
            end do
        else
           !+ Pass b0
            do i=1, ceiling(t/dt)
                call f(ax, ay, vx, vy, d, bc, dx, dy, kx, ky, ax0=ax0, ay0=ay0)
                ax = ax + dt*kx
                ay = ay + dt*ky
            end do
        end if
    end subroutine euler


    !+ Heun solver
    subroutine heun(ax, ay, vx, vy, d, bc, dx, dy, dt, t, ax0, ay0)
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy, dt, t
        real, intent(inout) :: ax(:,:), ay(:,:)
        integer, intent(in) :: bc
        real, optional, intent(in) :: ax0(:,:), ay0(:,:)
        real :: k1x(nx, ny+1), k1y(nx+1, ny), k2x(nx, ny+1), k2y(nx+1, ny)
        integer :: i

        if (present(ax0) .and. present(ay0)) then
            !+Pass b0
            do i=1, ceiling(t/dt)
                call f(ax, ay, vx, vy, d, bc, dx, dy, k1x, k1y, ax0=ax0, ay0=ay0)
                call f(ax + dt*k1x, ay + dt*k1y, vx, vy, d, bc, dx, dy, k2x, k2y, ax0=ax0, ay0=ay0)
                ax = ax + 0.5*dt*(k1x+k2x)
                ay = ay + 0.5*dt*(k1y+k2y)
            end do
        else
            !+ Don't pass b0
            do i=1, ceiling(t/dt)
                call f(ax, ay, vx, vy, d, bc, dx, dy, k1x, k1y)
                call f(ax + dt*k1x, ay + dt*k1y, vx, vy, d, bc, dx, dy, k2x, k2y)
                ax = ax + 0.5*dt*(k1x+k2x)
                ay = ay + 0.5*dt*(k1y+k2y)
            end do
        endif
    end subroutine heun
end module cft_kernel
