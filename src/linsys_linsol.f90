
submodule (linsys) linsys_linsol

  implicit none

contains

  !for Ax=b
  !x, b vectors
  module function linsol_vec(A,b) result(x)
    !Solves the linear system Ax=b for x.
    !---Inputs---
    !A: full rank square matrix (these are checked)
    !b: right hand side of linear system
    !---Outputs---
    !x: solution to linear system
    real, intent(in) :: A(:,:)
    real, intent(in) :: b(size(A,2))
    real, allocatable :: x(:)
    real, allocatable, dimension(:,:) :: L, U
    real, allocatable :: b_permuted(:), y(:)
    integer, allocatable :: p(:) !row pivot vector (such that PA=A(p,:))
    integer :: i, n

    if (.not. is_square(A,n)) then
       print *, "ERROR: cannot solve nonsquare linear system"
       return
    end if

    allocate(L(n,n), U(n,n), b_permuted(n), y(n), p(n))

    call lupp(A,L,U,p)

    if (minval(abs(get_diag(U))) == 0.0) then
       print *, "ERROR: linear system is singular"
       return
    end if
  
    b_permuted=b(p)
    y=fwd_sub(L,b_permuted) !solve L(Ux)=b_permuted for y=Ux
    x=bck_sub(U,y) !solve Ux=y=(Ux) for x
  end function linsol_vec


  !for AX=B
  !X, B matrices
  module function linsol_mat(A,B) result(X)
    !Solves the linear system AX=B for X.
    !---Inputs---
    !A: full rank square matrix (these are checked)
    !B: matrix right hand side of linear system
    !---Outputs---
    !X: matrix solution to "matrix linear system"
    real, intent(in) :: A(:,:)
    real, intent(in) :: B(:,:)
    real, allocatable :: X(:,:)
    real, allocatable, dimension(:,:) :: L, U
    real, allocatable :: b_cur(:), b_permuted(:), y(:)
    integer, allocatable :: p(:) !row pivot vectors (such that PA=A(p,:)
    integer :: i, n

    if (.not. is_square(A,n)) then
       print *, "ERROR: cannot solve nonsquare linear system"
       return
    end if

    if (size(B,1) /= n) then
       print *, "ERROR: dimension 1 of matrix A and matrix B do not match"
       return
    end if

    allocate(L(n,n), U(n,n), X(n,size(B,2)))
    allocate(b_cur(n), b_permuted(n), y(n), p(n))

    call lupp(A,L,U,p)

    if (minval(abs(get_diag(U))) == 0.0) then
       print *, "ERROR: linear system is singular"
       return
    end if
  
    do i=1,size(B,2)
       !pick out a column of B matrix, solve that linear system, insert into X
       b_cur=B(:,i) !pick out a column of B matrix
       b_permuted=b_cur(p)
       y=fwd_sub(L,b_permuted) !solve L(Ux)=b_permuted for y=Ux
       X(:,i)=bck_sub(U,y) !solve Ux=y=(Ux) for x
    end do
  end function linsol_mat

end submodule linsys_linsol
