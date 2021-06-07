
program test_linsys
  implicit none

  !Tests all routines in linsys.f90
  real :: abs_tol=1e-8
  real :: sor_tol=1e-8 !convergence tolerance for sor()
  integer :: max_iter=100 !maximum iterations for sor()

  print *
  print *, "---LINSYS TESTS---"
  call test_linsol(abs_tol)
  call test_sor(abs_tol,sor_tol,max_iter)
  print *
end program test_linsys


subroutine test_linsol(abs_tol)
  !tests linear system solver (LU version)
  use linsys

  implicit none
  real, intent(in) :: abs_tol
  real, allocatable :: A(:,:), b(:), x(:)
  real :: abs_err
  integer :: n=3

  !allocate and fill appropriately sized arrays
  allocate(A(n,n), b(n), x(n))
  call random_number(A)
  call random_number(b)

  !solve system
  x=linsol(A,b)

  abs_err=norm2(matmul(A,x)-b)
  if (abs_err > abs_tol) then
     print *, "  FAILED,    linsol (linear system solve)"
     print *, "    residual 2-norm: ", abs_err
  else
     print *, "  PASSED,    linsol (linear system solve)"
  end if
     
end subroutine test_linsol


subroutine test_sor(abs_tol,sor_tol,max_iter)
  !tests successive over-relaxation
  use linsys

  implicit none
  real, intent(in) :: abs_tol, sor_tol
  integer, intent(in) :: max_iter
  real, allocatable :: A(:,:), b(:), x(:), d(:)
  real :: abs_err
  integer :: n=3 !system size

  !allocate and fill appropriately sized arrays
  allocate(A(n,n), b(n), x(n), d(n))
  call random_number(A)
  call random_number(b)

  d=1.0
  A=A+set_diag(n*d) !ensure A is diagonally dominant by increasing diagonal

  x=sor(A,b,sor_tol,max_iter)

  abs_err=norm2(matmul(A,x)-b)
  if (abs_err > abs_tol) then
     print *, "  FAILED,    sor (successive over-relaxation)"
     print *, "    residual 2-norm: ", abs_err
  else
     print *, "  PASSED,    sor (successive over-relaxation)"
  end if

end subroutine test_sor


