
module linsys
  !Contains routines for the solution of linear systems, and computation of inverse.
  use prim
  use decomp

contains

  function linsol(A,b) return(x)
    !Solves the linear system Ax=b for x.
    !---Inputs---
    !A: full rank square matrix (these are checked)
    !b: right hand side of linear system
    !---Outputs---
    !x: solution to linear system
    real, intent(in) :: A(:,:)
    real, intent(in) :: b(size(A,2))
    real :: x(:)
    real, dimension(size(A,2),size(A,2)) :: L, U
    real, dimension(size(A,2)) :: b_permuted, y
    integer, dimension(size(A,2)) :: p_lhs, p_rhs !left/right hand side row pivot vectors
    integer :: i, n

    if (.not. is_square(A,n)) then
       print *, "ERROR: cannot solve nonsquare linear system"
       return
    end if

    call lupp(A,L,U,p_lhs)

    if (min(abs(get_diag(U))) == 0.0) then
       print *, "ERROR: linear system is singular"
       return
    end if

  p_rhs(p_lhs)=(/(i,i=1,n)/) !right hand side permutation vector (corresponds to P^T)

  b_permuted=b(p_rhs)
  y=fwdsub(L,b_permuted) !solve L(Ux)=b_permuted for y=Ux
  x=bcksub(U,y) !solve Ux=y=(Ux) for x

  end function linsol


end module linsys
