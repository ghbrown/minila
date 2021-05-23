
module prim
!primitive linear algebraic operations not included in standard Fortran
implicit none

contains

  subroutine dispmat(A)
    !display matrix
    !---Inputs---
    !A: a matrix
    !---Outputs---
    !NONE: just prints out matrix
    implicit none

    real, intent(in) :: A(:,:)
    integer :: i

    do i=1,size(A,2)
       print *, A(i,:)
    end do
  end subroutine dispmat
  

  function rop(u1,u2) result(u1u2T)
    !rank one product
    !---Inputs---
    !u1: order 1 array (vector)
    !u2: order 1 array (vector)
    !---Outputs---
    !u1u2T: order 2 array equal to (u1)transpose(u2)
    
    implicit none
    real, intent(in) :: u1(:)
    real, intent(in) :: u2(:)
    real :: u1u2T(size(u1),size(u2)) !variable for rank 1 product
    real :: u1m(size(u1),1) !variable for column matrix of u1
    real :: u2m(1,size(u2)) !variable for row matrix form of u2

    u1m=reshape(u1,(/size(u1),1/)) !shape first vector into column matrix
    u2m=reshape(u2,(/1,size(u2)/)) !shape second vector into row matrix

    u1u2T=matmul(u1m,u2m) !compute rank 1 product
  end function rop

end module prim


