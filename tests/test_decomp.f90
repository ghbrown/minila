
program test_decomp
  implicit none

  !had to depart from 1e-10 here for LUP
  real :: abs_tol=1e-8

  print *, "---DECOMPOSITION TESTS---"
  call test_lupp(abs_tol)

end program test_decomp


subroutine test_lupp(abs_tol)
  use prim
  use decomp

  implicit none
  real, intent(in) :: abs_tol
  !real, dimension(4,4):: A=transpose(reshape((/2.,0.,2.,0.6,3.,3.,4.,-2.,5.,5.,4.,2.,-1.,-2.,3.4,-1./),(/4,4/)))
  !the A initilized above is from Cormen's "Introduction to Algorithms"
  !and is used for testing when new features are added to lupp
  !the call to random_number below is used for general testing after
  !lupp has been verified on the Cormen example
  !real, dimension(4,4) :: L, U
  !integer, dimension(4) :: p, pT !permutation vectors corresponding to matrices P and P^T
  real, allocatable, dimension(:,:) :: A, L, U
  integer, allocatable, dimension(:) :: p, pT !permutation vectors corresponding to matrices P and P^T
  integer :: i, n=3
  real :: res_mean

  allocate(A(n,n), L(n,n), U(n,n))
  allocate(p(n),pT(n))

  call random_number(A)
  call lupp(A,L,U,p)

  pT(p)=(/(i,i=1,size(p,1))/) !permutation vector corresponding to P^T

  res_mean=sum(abs(A-matmul(L(pT,:),U)))/size(A) !mean of the absolute value of the matrix residual

  if (res_mean > abs_tol) then
     print *, "  FAILED,    lupp (LU, partial pivoting)"
     print *, "    mean element residual:    ", res_mean
  else
     print *, "  PASSED,    lupp (LU, partial pivoting)"
  end if
  
end subroutine test_lupp


