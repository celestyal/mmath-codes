program tests_find_dt
    use cft_kernel, only : find_dt

    implicit none

    ! Perform tests
    call test_find_dt()

    contains
        subroutine test_find_dt()
            logical :: passed, test_results(1)

            passed = .true.

            test_results(1) = test_standard()

            !+ If any test failed, exit with code 1
            if(any(test_results) .eqv. .false.) stop 1
        end subroutine test_find_dt

        function test_standard() result(test_result)
            real, allocatable :: vx(:,:), vy(:,:)
            real :: test_output
            real, parameter :: expected_output=30.0
            logical :: test_result
            character(128) :: message
            
            allocate(vx(100,100), vy(100, 100))
            vx(:,:) = 2.0
            vy(:,:) = 0.0
            test_output = find_dt(vx, vy, 300.0, 300.0, 0.0)
            test_result = test_output == expected_output

            if (test_result .eqv. .false.) then
                write(message, '(A,F6.2,A,F6.2)') "find_dt test Failed: Expected value of", expected_output, " but got", test_output
            else
                write(message, '(A)') 'Test passed'
            end if
            print *, message
        end function test_standard
end program tests_find_dt
