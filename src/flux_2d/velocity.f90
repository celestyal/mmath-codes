module velocity

    implicit none

contains


    subroutine apply_convergent_flow(velocity, speed, convergence_point, coordinate, nx)
        real, intent(inout) :: velocity(:,:)
        real, intent(in) :: speed, convergence_point, coordinate(:)
        integer, intent(in) :: nx
        integer :: i

        do concurrent (i = 1:nx+1, coordinate(i) .lt. convergence_point)
        velocity(:,i) = velocity(:,i) + speed
        end do
        do concurrent (i = 1:nx+1, coordinate(i) .gt. convergence_point)
        velocity(:,i) = velocity(:,i) - speed
        end do

    end subroutine apply_convergent_flow


    subroutine apply_differential_flow(velocity, omega, coordinate, nx)
        real, intent(inout) :: velocity(:,:)
        real, intent(in) :: omega, coordinate(:)
        integer, intent(in) :: nx
        integer :: i

        ! to finish implementing

        do concurrent (i=1:nx+1)
        velocity(:,i) = -omega*(coordinate(i)-0.0)
        end do
    end subroutine apply_differential_flow


    subroutine apply_meridional_flow(velocity, speed)
        real, intent(inout) :: velocity(:,:), speed

        ! to finish implementing

        velocity = velocity + speed
    end subroutine apply_meridional_flow


    subroutine apply_velocity_scaling(diffusion, x, y, length)
        real, intent(inout) :: diffusion, x(:,:), y(:,:), length

        diffusion = diffusion / (length**2)
        x = x / length
        y = y / length
    end subroutine apply_velocity_scaling

end module velocity
