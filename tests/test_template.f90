
program test_set
  implicit none

  !Tests all routines in a set by calling relevant routine
  !specific tests.
  real :: abs_tol=1e-10

  print *
  print *, "--- ___ TESTS---"
  !call test subroutines
  print *
end program test_


subroutine test_routine(abs_tol)
  !Tests a single function from the set
  use set

  implicit none
end subroutine test_routine

