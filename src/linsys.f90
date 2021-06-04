
module linsys
  !Contains routines for the solution of linear systems, and computation of inverse.
  use prim
  use decomp
  implicit none

contains

  function linsol(A,b) result(x)
    !Solves the linear system Ax=b for x.
    !---Inputs---
    !A: full rank square matrix (these are checked)
    !b: right hand side of linear system
    !---Outputs---
    !x: solution to linear system
    real, intent(in) :: A(:,:)
    real, intent(in) :: b(size(A,2))
    real :: x(size(A,2))
    real, dimension(size(A,2),size(A,2)) :: L, U
    real, dimension(size(A,2)) :: b_permuted, y
    integer, dimension(size(A,2)) :: p !row pivot vectors (such that PA=A(p,:)
    integer :: i, n

    if (.not. is_square(A,n)) then
       print *, "ERROR: cannot solve nonsquare linear system"
       return
    end if

    call lupp(A,L,U,p)

    if (minval(abs(get_diag(U))) == 0.0) then
       print *, "ERROR: linear system is singular"
       return
    end if
  
    b_permuted=b(p)
    y=fwd_sub(L,b_permuted) !solve L(Ux)=b_permuted for y=Ux
    x=bck_sub(U,y) !solve Ux=y=(Ux) for x

  end function linsol


  function inv(A) result(A_inv)
    !Computes the inverse of a square matrix
    implicit none
    real, intent(in) :: A(:,:)
    real, allocatable :: A_inv(:,:)
    integer :: n

    if (.not. is_square(A,n)) then
       print *, "ERROR: nonsquare matrices have no inverse"
       print *, "       (perhaps you'd like a psuedoinverse?)"
       return
    end if

    allocate(A_inv(n,n))

  end function inv


end module linsys
