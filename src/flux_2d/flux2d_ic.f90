module vector_potential
    implicit none
contains
    subroutine current_sheet(ax, ay, x, nx)
        integer :: i
        integer, intent(in) :: nx
        real, intent(out) :: ax(:,:), ay(:,:)
        real, intent(in) :: x(:)
        real :: b0

        ax = 0.0

        b0 = 10.0
        do concurrent (i=1:nx+1)
        ay(i,:) = log(cosh(b0 * x(i)))
        end do
    end subroutine current_sheet


    subroutine advection_test(ax, ay, x, y, nx, ny)
        integer :: i
        integer, intent(in) :: nx, ny
        real, intent(out) :: ax(:,:), ay(:,:)
        real, intent(in) :: x(:), y(:)
        real :: b0

        b0 = 10.0
        do concurrent (i=1:nx+1)
        ay(i,:) = 0.5*b0*x(i)
        end do
        do concurrent (i=1:ny+1)
        ax(:,i) = -0.5*b0*y(i)
        end do
    end subroutine advection_test


    subroutine bipole_no_tilt(ax, ay, x, y, nx, ny)
        real :: rho, b0
        integer :: i, j
        integer, intent(in) :: nx, ny
        real, intent(out) :: ax(:,:), ay(:,:)
        real, intent(in) :: x(:), y(:)
        real, dimension(nx+1, ny) :: epsilon

        ! peak field strength
        b0 = 10.0

        ! separation between peaks
        rho = 0.20

        ax(:,:) = 0.0

        do concurrent (i=1:nx+1, j=1:ny)
        epsilon(i,j) = (0.5 * (x(i)**2) + y(j)**2) / (rho**2)
        end do

        ay = b0*rho*exp(0.5-epsilon)
    end subroutine bipole_no_tilt
end module vector_potential
