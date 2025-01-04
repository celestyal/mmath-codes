module mag_grid
    implicit none
contains
    subroutine populate_grid_coordinates(axx, axy, ayx, ayy, gridx, gridy, &
            bx, by, xmin, ymin, dx, dy, nx, ny)
        real, intent(out) :: axx(:), axy(:), ayx(:), ayy(:), gridx(:), &
            gridy(:), bx(:), by(:)
        real, intent(in) :: xmin, ymin, dx, dy
        integer, intent(in) :: nx, ny

        gridx =  velocity_x_coordinate(xmin, dx, nx)
        gridy =  velocity_y_coordinate(ymin, dy, ny)
        bx =  b_x_coordinate(xmin, dx, nx)
        by = b_y_coordinate(ymin, dy, ny)
        axx =  ajx_x_coordinate(xmin, dx, nx)
        axy = ajx_y_coordinate(ymin, dy, ny)
        ayx = ajy_x_coordinate(xmin, dx, nx)
        ayy = ajy_y_coordinate(ymin, dy, ny)
    end subroutine populate_grid_coordinates


    pure function velocity_x_coordinate(xmin, dx, nx) result(x)
        real, intent(in) :: xmin, dx
        integer, intent(in) :: nx
        real :: x(nx+1)
        integer :: i

        do concurrent (i=1:nx+1)
        x(i) = xmin + (dx * (i-1))
        end do
    end function velocity_x_coordinate


    pure function velocity_y_coordinate(ymin, dy, ny) result(y)
        real, intent(in) :: ymin, dy
        integer, intent(in) :: ny
        real :: y(ny+1)
        integer :: i

        do concurrent (i=1:ny+1)
        y(i) = ymin + (dy * (i-1))
        end do
    end function velocity_y_coordinate


    pure function b_x_coordinate(xmin, dx, nx) result(x)
        real, intent(in) :: xmin, dx
        integer, intent(in) :: nx
        real :: x(nx+2)

        x(2:nx+2) = velocity_x_coordinate(xmin, dx, nx)
        x(1) = xmin - (dx/2)
        x(2:nx+2) = x(2:nx+2) + 0.5*dx
    end function b_x_coordinate


    pure function b_y_coordinate(ymin, dy, ny) result(y)
        real, intent(in) :: ymin, dy
        integer, intent(in) :: ny
        real :: y(ny+2)

        y(2:ny+2) = velocity_y_coordinate(ymin, dy, ny)
        y(1) = ymin - dy/2
        y(2:ny+2) = y(2:ny+2) + 0.5*dy
    end function b_y_coordinate


    pure function ajx_x_coordinate(xmin, dx, nx) result(x)
        real, intent(in) :: xmin, dx
        integer, intent(in) :: nx
        real:: x(nx), x_temp(nx+1)

        x_temp = velocity_x_coordinate(xmin, dx, nx)
        x = x_temp(1:nx) + 0.5*dx
    end function ajx_x_coordinate


    pure function ajx_y_coordinate(ymin, dy, ny) result(y)
        real, intent(in) :: ymin, dy
        integer, intent(in) :: ny
        real :: y(ny+1)

        y =  velocity_y_coordinate(ymin, dy, ny)
    end function ajx_y_coordinate


    pure function ajy_x_coordinate(xmin, dx, nx) result(x)
        real, intent(in) :: xmin, dx
        integer, intent(in) :: nx
        real:: x(nx+1)

        x = velocity_x_coordinate(xmin, dx, nx)
    end function ajy_x_coordinate


    pure function ajy_y_coordinate(ymin, dy, ny) result(y)
        real, intent(in) :: ymin, dy
        integer, intent(in) :: ny
        real :: y(ny), y_temp(ny+1)

        y_temp = velocity_y_coordinate(ymin, dy, ny)
        y = y_temp(1:ny) + 0.5*dy
    end function ajy_y_coordinate
end module mag_grid
