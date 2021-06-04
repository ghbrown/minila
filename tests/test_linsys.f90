
program test_linsys
  implicit none

  !Tests all routines in linsys.f90
  real :: abs_tol=1e-8

  print *
  print *, "---LINSYS TESTS---"
  call test_linsol(abs_tol)
  print *
end program test_linsys


subroutine test_linsol(abs_tol)
  !Tests a single function from the set
  use linsys

  implicit none
  real, intent(in) :: abs_tol
  real, allocatable :: A(:,:)
  real, allocatable :: b(:)
  real, allocatable :: x(:)
  integer :: n=3

  !allocate and set up appropriately sized arrays
  allocate(A(n,n), b(n), x(n))

  call random_number(A)
  call random_number(b)

  !solve system
  x=linsol(A,b)

  if (norm2(matmul(A,x)-b) > abs_tol) then
     print *, "  FAILED,    linsol (linear system solve)"
     print *, "    residual: ", norm2(matmul(A,x)-b)
  else
     print *, "  PASSED,    linsol (linear system solve)"
  end if
     
end subroutine test_linsol


