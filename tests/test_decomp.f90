
program test_decomp
  implicit none

  real :: abs_tol=1e-10

  print *, "---DECOMPOSITION TESTS---"
  call test_lupp(abs_tol)

end program test_decomp


subroutine test_lupp(abs_tol)
  use prim
  use decomp

  implicit none
  real, intent(in) :: abs_tol
  real, dimension(5,5) :: A, L, U
  integer, dimension(5) :: P
  real :: res_mean
  call random_number(A)
  call lupp(A,L,U,P)
  !need to make matrix out of P and premultiply with L
  res_mean=sum(abs(A-matmul(L,U)))/size(A)

  if (res_mean > abs_tol) then
     print *, "  FAILED,    lupp (LU, partial pivoting)"
     print *, res_mean
  else
     print *, "  PASSED,    lupp (LU, partial pivoting)"
  end if
  
end subroutine test_lupp


