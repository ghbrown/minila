
module linsys
  !Contains routines for the solution of linear systems, and computation of inverse.
  use prim
  use decomp
  implicit none

contains

  !DETERMINISTIC METHODS
  !LU decomposition based linear solver
  interface
     !for vector unknown x
     function linsol(A,b) result(x)
        real, intent(in) :: A(:,:)
        real, intent(in) :: b(size(A,2))
        real, allocatable :: x(:)
     end function linsol

     !for matrix unknown X
     function linsol(A,B) result(X)
        real, intent(in) :: A(:,:)
        real, intent(in) :: B(:,:)
        real, allocatable :: X(:,:)
     end function linsol
  end interface

  !inverse
  interface
     function inv(A) result(A_inv)
       real, intent(in) :: A(:,:)
       real, allocatable :: A_inv(:,:)
     end function inv
  end interface


  !ITERATIVE METHODS
  !successive over-relaxation
  interface
     function sor(A,b,rel_tol,max_iter) result(x)
        real, intent(in) :: A(:,:), rel_tol
        real, intent(in) :: b(:)
        integer, intent(in) :: max_iter
        real, allocatable :: x(:)
     end function sor
  end interface
   
end module linsys
