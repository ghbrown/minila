
module linsys
  !routines for the solution of linear systems, and computation of inverse
  use prim
  use decomp
  implicit none

  private
  public :: linsol
  public :: inverse
  public :: sor

  !DETERMINISTIC METHODS
  !LU decomposition based linear solver
  interface linsol
     !for vector unknown x
     module function linsol_vec(A,b) result(x)
        real, intent(in) :: A(:,:)
        real, intent(in) :: b(size(A,2))
        real, allocatable :: x(:)
     end function linsol_vec

     !for matrix unknown X
     module function linsol_mat(A,B) result(X)
        real, intent(in) :: A(:,:)
        real, intent(in) :: B(:,:)
        real, allocatable :: X(:,:)
      end function linsol_mat
  end interface linsol

  !inverse
  interface inverse
     module function inverse(A) result(A_inv)
       real, intent(in) :: A(:,:)
       real, allocatable :: A_inv(:,:)
     end function inverse
  end interface inverse


  !ITERATIVE METHODS
  !successive over-relaxation
  interface sor
     module function sor(A,b,rel_tol,max_iter) result(x)
        real, intent(in) :: A(:,:), b(:), rel_tol
        integer, intent(in) :: max_iter
        real, allocatable :: x(:)
     end function sor
  end interface sor
   
end module linsys
