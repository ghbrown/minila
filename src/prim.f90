
module prim
  !primitive linear algebraic operations not included in standard Fortran
  !disp_mat(A): displays matrix A
  !strang(T,b): solves triangular system Tx=b
  implicit none

contains

  subroutine disp_mat(A)
    !display matrix
    !---Inputs---
    !A: matrix
    !---Outputs---
    !NONE: just prints out matrix
    implicit none
    real, intent(in) :: A(:,:)
    integer :: i

    do i=1,size(A,2)
       print *, A(i,:)
    end do
  end subroutine disp_mat


  function is_square(A,n) result(tof)
    !Determines if a matrix is square, and if so what it's
    !dimension is.
    !---Inputs---
    !A: a matrix (to be tested for triangularity)
    !---Outputs---
    !tof: "true or false" logical on whether A is square
    !n: dim(A) if A is square
    implicit none
    real, intent(in) :: A(:,:)
    integer, intent(out) :: n
    logical :: tof

    if (size(A,1) /= size(A,2)) then
       tof=.false. !T is not square
    else
       tof=.true.
       n=size(A,2) !define system size
    end if
  end function is_square

  
  function is_triangular(A,type) result(tof)
    !tests matrix for triangularity and outputs character
    !corresponding to triangularity type
    !---Inputs---
    !A: a matrix (to be tested for triangularity)
    !---Outputs---
    !tof: "true or false" logical on whether A is triangular
    !type: character defining type of triangularity:
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

    if (.not. is_square(A,n)) then
       tof=.false. !T is not square
       type="q"
       return
    else
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
    end if

  end function is_triangular


  function get_diag(A) result(diag_vec)
    !Returns diagonal of a square matrix.
    !---Inputs---
    !A: matrix
    !---Outputs---
    !diag_vec: diagonal of matrix A
    implicit none
    real, intent(in) :: A(:,:)
    real :: diag_vec(size(A,2))
    integer :: n, i

    if (.not. is_square(A,n)) then
       print *, "ERROR: non-square matrices do not have a sensible diagonal"
    else
       do i=1,n
          diag_vec(i)=A(i,i)
       end do
    end if
  end function get_diag


  function set_diag(diag_vec) result(A)
    !generates square matrix with specified diagonal
    !---Inputs---
    !diag_vec: vector of diagonal elements
    !---Outputs---
    !A: matrix with diagonal diag_vec
    implicit none
    real, intent(in) :: diag_vec(:)
    real :: A(size(diag_vec),size(diag_vec))
    integer :: n, i
    n=size(diag_vec)
    A=0.0
    do i=1,n
        A(i,i)=diag_vec(i)
    end do
  end function set_diag


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


  function diag_sub(D,b) result(x)
    !solve diagonal system Dx=b, solve using specific O(n) method
    !NOTE: does not check if D is diagonal (does check squareness)
    !---Inputs---
    !D: diagonal matrix
    !b: right hand side vector
    !---Outputs---
    !x: solution to diagonal system
    implicit none
    real, intent(in) :: D(:,:)
    real, intent(in) :: b(size(D,2))
    real :: x(size(D,2))
    integer :: n !system size/dimension
    integer :: i
    !LOOK FOR VECTORIZED VERSION LIKE:
    !x=b/diag(T)
    if (.not. is_square(D,n)) then
       print *, "ERROR: solving diagonal system requires square matrix"
    else
       x=b/get_diag(D)
    end if
  end function diag_sub


  function bck_sub(T,b) result(x)
    !Apply back substitution to solve linear system Tx=b
    !with upper triangular T
    !NOTE: does not check triangularity of T (does check squareness)
    !---Inputs---
    !T: upper triangular matrix
    !b: right hand side vector
    !---Outputs---
    !x: solution to triangular system
    implicit none
    real, intent(in) :: T(:,:)
    real, intent(in) :: b(size(T,2))
    real :: x(size(T,2))
    integer :: n !system size/dimension
    integer :: i
    if (.not. is_square(T,n)) then
       print *, "ERROR: back substitution requires square matrix"
    else
       do i=n,1,-1
           x(i)=(1 / T(i,i)) * ( b(i) - dot_product(T(i,i+1:),x(i+1:)) )
       end do
    end if
  end function bck_sub


  function fwd_sub(T,b) result(x)
    !Apply forward substitution to solve linear system Tx=b
    !with lower triangular T
    !NOTE: does not check triangularity of T (does check squareness)
    !---Inputs---
    !T: lower triangular matrix
    !b: right hand side vector
    !---Outputs---
    !x: solution to triangular system
    implicit none
    real, intent(in) :: T(:,:)
    real, intent(in) :: b(size(T,2))
    real :: x(size(T,2))
    integer :: n !system size/dimension
    integer :: i
    if (.not. is_square(T,n)) then
       print *, "ERROR: forward substitution requires square matrix"
    else
       do i=1,n
          x(i)=(1 / T(i,i)) * ( b(i) - dot_product(T(i,1:i-1),x(1:i-1)) )
       end do
    end if
  end function fwd_sub


  function strang(T,b) result(x)
    !solve triangular system, detects type of triangularity
    !convenient wrapper and safeguarder for (diag/bck/fwd)sub
    !NOTE: does an O(n^2) check for type of triangularity,
    !      if performance is critical use the *sub functions
    !      directly
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
       x=diag_sub(T,b) !use situational O(n) method to solve.
    else if (type == "l") then
       !T is lower triangular
       x=bck_sub(T,b) !use back substitution
    else
       !T is triangular, but not lower triangular (must be upper)
       x=fwd_sub(T,b) !use forward substitution
    end if

  end function strang
  

end module prim


