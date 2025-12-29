program tests_find_dt
    use cft_kernel, only : find_dt

    implicit none

    ! Test result variables
    integer, parameter :: num_tests = 1
    logical :: passed(num_tests)

    ! Perform tests
    passed(1) = test_find_dt()

    !+ If any test failed, exit with code 1
    if(any(passed) .eqv. .false.) stop 1

    contains
        function test_find_dt() result(passed)
            logical :: passed, test_result(1)
            real :: test_output
            real, allocatable :: vx(:,:), vy(:,:)
            character(128) :: message

            passed = .true.

            !+ First test
            allocate(vx(100,100), vy(100, 100))
            vx(:,:) = 2.0
            vy(:,:) = 0.0
            test_output = find_dt(vx, vy, 300.0, 300.0, 0.0)
            test_result(1) = test_output == 30.0
            if (test_result(1) .eqv. .true.) then
                write(message, '(A)') 'Test passed'
            else
                write(message, '(A,F6.2,A,F6.2)') "find_dt test Failed: Expected value of ", 30.0, " but got", test_output
            end if
            print *, message
            
            !+ Check if any test failed
            if(any(test_result) .eqv. .false.) passed = .false.
        end function test_find_dt
end program tests_find_dt
