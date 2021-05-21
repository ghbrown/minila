
module primitives
!primitive linear algebraic operations not included in standard Fortran
implicit none

contains

  subroutine r1p(u1,u2,u1Tu2)
    !rank one product
    !---Inputs---
    !u1: order 1 array (vector)
    !u2: order 1 array (vector)
    !---Outputs---
    !u1Tu2: order 2 array equal to transpose(u1)u2
    
    implicit none
    real, intent(in) :: u1(:)
    real, intent(in) :: u2(:)
    real :: u1m(size(u1),1) !variable for column matrix of u1
    real :: u2m(1,size(u2)) !variable for row matrix form of u2
    real, intent(out) :: u1Tu2(size(u1),size(u2)) !variable for rank 1 product

    u1m=reshape(u1,(/size(u1),1/)) !shape first vector into column matrix
    u2m=reshape(u2,(/1,size(u2)/)) !shape second vector into row matrix

    u1Tu2=matmul(u1m,u2m) !compute rank 1 product
  end subroutine r1p

end module primitives


