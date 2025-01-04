module numerical_solver
    use file_io
    implicit none

    ! quantitys of interest, globally accessible
    real, allocatable :: b(:,:), jx(:,:), jy(:,:)

    ! these are temporary arrays used for evolution, since we want to avoid updating the
    ! actual quantities.
    real, private, allocatable, dimension(:,:) :: b_corner, vxx, vxy, vxx_yrib, &
        vxy_xrib
contains

    subroutine init_farrays(nx, ny)
        integer, intent(in) :: nx, ny

        ! initialise the arrays outside of the loop so that memory isn't being allocated and deallocated during loop.
        allocate(b_corner(nx+1, ny+1))
        allocate(vxx(nx+1, ny+1))
        allocate(vxy(nx+1, ny+1))
        allocate(vxx_yrib(nx+1, ny))
        allocate(vxy_xrib(nx, ny+1))
        allocate(b(nx+2,ny+2))
        allocate(jx(nx,ny+1))
        allocate(jy(nx+1,ny))
    end subroutine init_farrays


    subroutine f(ax, ay, vx, vy, d, bc, dx, dy, nx, ny, dax, day)
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy
        integer, intent(in) :: nx, ny, bc
        real, intent(in) :: ax(:,:), ay(:,:)
        real, intent(out) ::  dax(:,:), day(:,:)

        call calculate_b(ax, ay, bc, dx, dy, nx, ny)

        ! 4 point average for b at corner
        b_corner = (b(1:nx+1, 1:ny+1) + b(1:nx+1, 2:ny+2) + b(2:nx+2, 1:ny+1) + b(2:nx+2, 2:ny+2)) / 4.0

        ! on the corners
        vxx = b_corner * vx
        vxy = b_corner * vy

        ! on the locations of ax and ay
        vxx_yrib = (vxx(1:nx+1, 1:ny) + vxx(1:nx+1, 2:ny+1)) / 2.0
        vxy_xrib = (vxy(1:nx, 1:ny+1) + vxy(2:nx+1, 1:ny+1)) / 2.0

        ! partial derivatives of jx and jy
        call calculate_j(dx, dy, nx, ny)

        ! induction equation
        dax = vxy_xrib - d*jx
        day = -vxx_yrib - d*jy
    end subroutine f


    subroutine calculate_b(ax, ay, bc, dx, dy, nx, ny)
        real, intent(in) :: ax(:,:), ay(:,:), dx, dy
        integer, intent(in) :: bc, nx, ny

        b(2:nx+1,2:ny+1) = (ay(2:nx+1, 1:ny)-ay(1:nx, 1:ny)) / dx - (ax(1:nx, 2:ny+1)-ax(1:nx, 1:ny)) / dy
        if (bc .ne. 0) call apply_bc(bc, nx, ny)
    end subroutine calculate_b


    subroutine apply_bc(bc, nx, ny)
        integer, intent(in) :: bc, nx, ny
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
    end subroutine apply_bc


    subroutine calculate_j(dx, dy, nx, ny)
        real, intent(in) :: dx, dy
        integer, intent(in) :: nx, ny

        jx = (b(2:nx+1, 2:ny+2) - b(2:nx+1, 1:ny+1)) / dy
        jy =-(b(2:nx+2, 2:ny+1) - b(1:nx+1, 2:ny+1)) / dx
    end subroutine calculate_j


    subroutine rk4(ax, ay, vx, vy, d, bc, dx, dy, dt, n_steps, nx, ny, start, n_saves, id)
        real, intent(inout) :: ax(:,:), ay(:,:)
        integer, intent(in) :: bc, nx, ny, n_steps, n_saves
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy, dt, start
        character(512), intent(in) :: id
        real, dimension(nx, ny+1) :: x1, x2, x3, x4
        real, dimension(nx+1, ny) :: y1, y2, y3, y4
        integer :: i

        do i=1, n_steps
            call f(ax, ay, vx, vy, d, bc, dx, dy, nx, ny, x1, y1)
            call f(ax + 0.5*dt*x1, ay + 0.5*dt*y1, vx, vy, d, bc, dx, dy, nx, ny, x2, y2)
            call f(ax + 0.5*dt*x2, ay + 0.5*dt*y2, vx, vy, d, bc, dx, dy, nx, ny, x3, y3)
            call f(ax + dt*x3, ay + dt*y3, vx, vy, d, bc, dx, dy, nx, ny, x4, y4)
            ax = ax + (dt/6)*(x1 + 2*x2 + 2*x3 + x4)
            ay = ay + (dt/6)*(y1 + 2*y2 + 2*y3 + y4)

            call solver_progress(i, n_steps, start, n_saves, id, ax, ay, dx, dy, nx, ny, dt)
        end do
    end subroutine rk4


    subroutine euler(ax, ay, vx, vy, d, bc, dx, dy, dt, n_steps, nx, ny, start, n_saves, id)
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy, dt, start
        real, intent(inout) :: ax(:,:), ay(:,:)
        integer, intent(in) :: bc, nx, ny, n_steps, n_saves
        character(512), intent(in) :: id
        real :: kx(nx, ny+1), ky(nx+1, ny)
        integer :: i

        do i=1, n_steps
            call f(ax, ay, vx, vy, d, bc, dx, dy, nx, ny, kx, ky)
            ax = ax + dt*kx
            ay = ay + dt*ky

            call solver_progress(i, n_steps, start, n_saves, id, ax, ay, dx, dy, nx, ny, dt)
        end do
    end subroutine euler


    subroutine heun(ax, ay, vx, vy, d, bc, dx, dy, dt, n_steps, nx, ny, start, n_saves, id)
        real, intent(in) :: vx(:,:), vy(:,:), d, dx, dy, dt, start
        real, intent(inout) :: ax(:,:), ay(:,:)
        integer, intent(in) :: bc, nx, ny, n_steps, n_saves
        character(512), intent(in) :: id
        real :: k1x(nx, ny+1), k1y(nx+1, ny), k2x(nx, ny+1), k2y(nx+1, ny)
        integer :: i

        do i=1, n_steps
            call f(ax, ay, vx, vy, d, bc, dx, dy, nx, ny, k1x, k1y)
            call f(ax + dt*k1x, ay + dt*k1y, vx, vy, d, bc, dx, dy, nx, ny, k2x, k2y)
            ax = ax + 0.5*dt*(k1x+k2x)
            ay = ay + 0.5*dt*(k1y+k2y)

            call solver_progress(i, n_steps, start, n_saves, id, ax, ay, dx, dy, nx, ny, dt)
        end do
    end subroutine heun


    subroutine solver_progress(i, n_steps, start, n_saves, id, ax, ay, dx, dy, nx, ny, dt)
        integer, intent(in) :: i, n_steps, n_saves, nx, ny
        real, intent(in) :: ax(:,:), ay(:,:), start, dx, dy, dt
        character(512), intent(in) :: id
        real :: progress
        real, allocatable :: save_points(:)
        integer :: k

        progress = (real(i)/real(n_steps))*100
        if (mod(i, ceiling(real(n_steps)/100.0)) == 0.0 ) then
            print *, "progress: ", nint(progress), "%"
        end if

        save_points = [(floor(real(n_steps)/(real(n_saves)-1.0))*k,k=1,n_saves-1)]

        if (n_saves .gt. 2) then
            if ((i == n_steps) .or. (any(save_points==i))) then
                ! save system
                call save_snapshot(.false., id, ax, ay, dx, dy, nx, ny)
                ! save corresponding timestamp
                call save_timestamp(id, (start+i*dt))
            end if
        end if
    end subroutine solver_progress


    subroutine save_snapshot(is_initial, id, ax, ay, dx, dy, nx, ny, vx, vy)
        real, intent(in) :: ax(:,:), ay(:,:), dx, dy
        real, intent(in), optional :: vx(:,:), vy(:,:)
        logical, intent(in) :: is_initial
        character(512), intent(in) :: id
        integer :: nx, ny
        character(6) :: mode

        if (is_initial) then
            mode = "asis"
            call export_data(id, vx, "vx.dat", trim(mode))
            call export_data(id, vy, "vy.dat", trim(mode))
        else
            mode = "append"
        end if

        call calculate_b(ax, ay, 0, dx, dy, nx, ny)
        call calculate_j(dx, dy, nx, ny)
        call export_data(id, b(2:nx+1,2:ny+1), "bz.dat", trim(mode)) ! exclude extrapolated points
        call export_data(id, ax, "ax.dat", trim(mode))
        call export_data(id, ay, "ay.dat", trim(mode))
        if (present(vx)) call export_data(id, vx, "vx.dat", trim(mode))
        if (present(vy)) call export_data(id, vy, "vy.dat", trim(mode))
        call export_data(id, jx, "jx.dat", trim(mode))
        call export_data(id, jy, "jy.dat", trim(mode))
    end subroutine save_snapshot
end module numerical_solver
