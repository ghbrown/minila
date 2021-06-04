module decomp
  !Contains subroutines for matrix decompositions: LU 
  !and eventually: QR, eigenvalue
  use prim
  implicit none

contains

  !TO DO:
  !fix LUPP decomposition to permit rectangular matrices
  !implement packed storage (must do AFTER the above)

  subroutine lupp(A,L,U,p)
    !Performs LU decomposition of A with partial pivoting.
    !---Inputs---
    !A: square matrix (this is checked)
    !---Outputs---
    !L: lower triangular matrix
    !U: upper triangular matrix
    !p: "left" pivot vector
    !NOTE: the vector p is returned such that P_mat A = A(p,:) = L U,
    !      where P_mat is the matrix constructable from P
    !      ****if p is to be used on the right hand side (to permute L,
    !      for example) one must compute the vector corresponding to
    !      transpose(P_mat):
    !         pT(p)=(/(i,i=1,size(p,1))/)
    !      which satisfies A = L(pT,:) U ***
    implicit none
    real, intent(in) :: A(:,:)
    real, intent(out), allocatable :: L(:,:), U(:,:)
    integer, intent(out), allocatable :: p(:) !vector used to keep track of column swaps
    real, allocatable :: S(:,:)
    integer :: i, n, cur_iter, max_row_rel

    !THIS ERROR IS ITSELF ERRONEOUS, LU DOES EXIST FOR NONSQUARE MATRICES
    !THIS MAY BE A PAIN TO IMPLEMENT (INDICES BELOW, ETC.)
    if (.not. is_square(A,n)) then !check if A is square
       print *, "ERROR: LU has not yet been implemented for"
       print *, "         nonsquare matrices"
       return !LU not implemented nonsquare matrices, exit early
    end if

    allocate(L(n,n), U(n,n), p(n), S(n,n))

    S=A
    L=0.0
    U=0.0
    p=(/(i,i=1,size(A,1))/)

    !rewrite using packed storage (the method described below)?
    !there may be a way to structure this such that it overwrites S into (L_strict + U, L_strict: strictly lower triangular, U: upper triangular with nontrivial diagonal elements) which would save memory (until somebody wants L and U)
    !this could be done by overwriting the non-sub-matrix part of S (one row-column shell per iteration)

    do cur_iter=1,n
      max_row_rel=maxloc(abs(S(cur_iter:,cur_iter)),1)-1 !get index of row with maximum leading entry (relative to cur_iter row)
      p((/cur_iter,cur_iter+max_row_rel/))=p((/cur_iter+max_row_rel,cur_iter/)) !perform swap of P elements in place
      S((/cur_iter+max_row_rel,cur_iter/),cur_iter:)=S((/cur_iter,cur_iter+max_row_rel/),cur_iter:) !perform swap of submatrix rows in place
      if ( S(cur_iter,cur_iter) == 0.0) then
         print *, "ERROR: Zero pivot, attempted to LU decompose singular matrix"
         return
      end if
      L((/cur_iter+max_row_rel,cur_iter/),:cur_iter)=L((/cur_iter,cur_iter+max_row_rel/),:cur_iter) !perform swap of rows in L which have already been computed 

      L(cur_iter,cur_iter)=1 !set diagonal element of L
      L(cur_iter+1:,cur_iter)=S(cur_iter+1:,cur_iter)/S(cur_iter,cur_iter)
      U(cur_iter,cur_iter:)=S(cur_iter,cur_iter:)
      S(cur_iter+1:,cur_iter+1:)=S(cur_iter+1:,cur_iter+1:)-rop(L(cur_iter+1:,cur_iter),U(cur_iter,cur_iter+1:)) !update Schur complement submatrix in S
    end do

  end subroutine lupp


end module decomp
