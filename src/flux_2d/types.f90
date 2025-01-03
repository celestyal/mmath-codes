module types
    implicit none

    type plane_component
        real, allocatable, dimension(:,:) :: x, y
    end type plane_component

    type coordinates
        real, allocatable, dimension(:) :: x, y
    end type coordinates

    type profile
        integer :: boundary_condition, potential_configuration, steps
        real :: convergence_speed, convergence_point, diffusivity, length_scale, &
            meridional_speed, differential_speed, simulation_end_time, simulation_start_time
    end type profile

    type debug
        logical :: enable_test, new_simulation
        integer :: numerical_solver, test_number, number_of_saves
    end type debug

    type grid
        real :: xmax, xmin, ymax, ymin
        integer :: nx, ny
    end type grid

end module types
