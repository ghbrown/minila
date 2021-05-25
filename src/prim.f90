
module prim
  !primitive linear algebraic operations not included in standard Fortran
  !dispmat(A): displays matrix A
  !strang(T,b): solves triangular system Tx=b
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


  !function should determine if it is upper or lower triangular and transpose lower triangular matrices (or use transposed indices) to always solve upper triangular system.
  function strang(T,b) result(x)
    !solve triangular system
    !---Inputs---
    !T: (upper OR lower) triangular matrix
    !b: right hand side vector
    !---Outputs---
    !x: solution to triangular system
    implicit none
    real, intent(in) :: T(:,:)
    real, intent(in) :: b(size(T,2))
    real :: x(size(T,2))
    integer :: n !system size/dimension
    integer :: i_forward !increasing index
    integer :: i !decreasing index computed from i_forward
    character :: type !will be set to "upper" or "lower"
    if (size(T,1) /= size(T,2)) then
       print *, "ERROR: Attempted to solve non-square triangular system."
       return
    else
       n=size(T,2)
    end if

    !Back substitution
    do i=n,1,-1
       x(i)=(1 / T(i,i)) * ( b(i) - dot_product(T(i,i+1:),x(i+1:)) )
    end do

  end function strang

  
  

end module prim


