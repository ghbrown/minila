
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

  function is_triangular(A,type) result(tof)
    !display matrix
    !---Inputs---
    !A: a matrix (to be tested for triangularity)
    !---Outputs---
    !type: character defining what type of triangular it is
    !      "q", not square
    !      "n", not triangular
    !      "u", upper triangular
    !      "l", lower triangular
    !      "d", diagonal (upper and lower triangular)
    real, intent(in) :: A(:,:)
    character, intent(out) :: type
    logical :: tof
    real :: sum_lower, sum_upper
    integer :: i, n

    if (size(A,1) /= size(A,2)) then
       tof=.false. !T is not square
       type="q"
    else
       n=size(A,2) !define system size
    end if

    sum_lower=0
    sum_upper=0
    do i=1,n
       sum_lower=sum_lower+sum(abs(A(i+1:,i)))
       sum_upper=sum_upper+sum(abs(A(1:i-1,i)))
    end do

    if ((sum_upper > 0) .and. (sum_lower > 0)) then
       tof=.false. !T is not triangular
       type="n"
    else if (sum_upper+sum_lower == 0.0) then
       type="d" !T is upper and lower triangular (diagonal)
       tof=.true. 
    else if (sum_upper > 0) then
       tof=.true. !sum of upper triangle is zero
       type="l" !T must be lower triangular
    else if (sum_lower > 0) then
       tof=.true. !sum of lower triangle is zero
       type="u" !T must be upper triangular
    end if

  end function is_triangular


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


  function strang(T,b) result(x)
    !solve triangular system, detects type of triangularity
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
    integer :: i
    character :: type !will be used as output of is_triangular

    if (.not. is_triangular(T,type)) then
       print *, "ERROR: Attempted to solve non-triangular system."
       return
    else
       n=size(T,2)
    end if

    if (type=="d") then
       !T is diagonal
       !Use O(n) method to solve.
       do i=1,n
          x(i)=b(i)/T(i,i) !
       end do
    else if (type == "l") then
       !T is lower triangular
       !use back substitution
       do i=n,1,-1
          x(i)=(1 / T(i,i)) * ( b(i) - dot_product(T(i,i+1:),x(i+1:)) )
       end do
    else
       !T is triangular, but not lower triangular (must be upper)
       !use forward substitution
       do i=1,n
          x(i)=(1 / T(i,i)) * ( b(i) - dot_product(T(i,1:i-1),x(1:i-1)) )
       end do
    end if

  end function strang
  

end module prim


