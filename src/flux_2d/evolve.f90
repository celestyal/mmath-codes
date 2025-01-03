module evolve
    implicit none

    ! these are temporary arrays used for evolution, since we want to avoid updating the
    ! actual quantities.
    real, private, allocatable, dimension(:,:) :: b_corner, vxx, vxy, vxx_yrib, &
        vxy_xrib

    subroutine init_farrays(nx, ny)
        integer, intent(in) :: nx, ny

        ! initialise the arrays outside of the loop so that memory isn't being allocated and deallocated during loop.
        allocate(b_corner(nx+1, ny+1))
        allocate(vxx(nx+1, ny+1))
        allocate(vxy(nx+1, ny+1))
        allocate(vxx_yrib(nx+1, ny))
        allocate(vxy_xrib(nx, ny+1))

    end subroutine init_farrays


    subroutine f(ax, ay, vx, vy, diffusivity, boundary_condition, dx, dy, nx, ny, dax, day)
        real, intent(in) :: vx(:,:), vy(:,:), diffusivity, dx, dy
        integer, intent(in) :: nx, ny, boundary_condition
        real, intent(in) :: ax(:,:), ay(:,:)
        real, intent(out) ::  dax(:,:), day(:,:)

        call calculate_b(ax, ay, boundary_condition, dx, dy, nx, ny)

        ! 4 point average for b at corner
        b_corner = ( b(1:nx+1, 1:ny+1) + b(1:nx+1, 2:ny+2) + b(2:nx+2, 1:ny+1) + b(2:nx+2, 2:ny+2) ) / 4.0

        ! on the corners
        vxx = b_corner * vx
        vxy = b_corner * vy

        ! on the locations of ax and ay
        vxx_yrib = ( vxx(1:nx+1, 1:ny) + vxx(1:nx+1, 2:ny+1) ) / 2.0
        vxy_xrib = ( vxy(1:nx, 1:ny+1) + vxy(2:nx+1, 1:ny+1) ) / 2.0

        ! partial derivatives of jx and jy
        call calculate_j(dx, dy, nx, ny)

        ! induction equation
        dax = vxy_xrib - diffusivity*jx
        day = -vxx_yrib - diffusivity*jy
    end subroutine f


    subroutine calculate_b(ax, ay, boundary_condition, dx, dy, nx, ny)
        real, intent(in) :: ax(:,:), ay(:,:), dx, dy
        integer, intent(in) :: boundary_condition, nx, ny

        b(2:nx+1,2:ny+1) = (ay(2:nx+1, 1:ny)-ay(1:nx, 1:ny)) / dx - (ax(1:nx, 2:ny+1)-ax(1:nx, 1:ny)) / dy

        call apply_boundary_condition(boundary_condition, nx, ny)

    end subroutine calculate_b


    subroutine apply_boundary_condition(boundary_condition, nx, ny)
        integer, intent(in) :: boundary_condition, nx, ny
        select case (boundary_condition)
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
    end subroutine apply_boundary_condition


    subroutine calculate_j(dx, dy, nx, ny)
        real, intent(in) :: dx, dy
        integer, intent(in) :: nx, ny

        jx = ( b(2:nx+1, 2:ny+2) - b(2:nx+1, 1:ny+1) ) / dy
        jy =-( b(2:nx+2, 2:ny+1) - b(1:nx+1, 2:ny+1) ) / dx
    end subroutine calculate_j
end module evolve
