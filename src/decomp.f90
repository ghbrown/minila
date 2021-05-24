module decomp
  !Contains subroutines for matrix decompositions: LU 
  !and eventually: QR, eigenvalue
  use prim
  implicit none

contains

  subroutine LU(A,L,U)
    implicit none
    real, intent(in) :: A(:,:)
    real, intent(out) :: L(size(A,1),size(A,1)), U(size(A,1),size(A,1))
    real :: S(size(A,1),size(A,1))
    integer :: n, cur_iter
    if (size(A,1) /= size(A,2)) then !check if A is square
       print *, 'ERROR: LU is invalid for nonsquare matrices'
       return !LU is not valid for nonsquare matrices, exit early
    end if

    n=size(A,1) !dimension of square system

    S=A
    L=0.0
    U=0.0
    !in principle it should be possible to do this whole algorithm without recursion by maintining an "A-sized" S throughout, but only using progressively less and less of it (shrinking the used submatrix) by also indexing it with cur_iter
    !this will take more memory and about twice as long, but it's infinitely more readable

    !there may also be a way to structure this such that it overwrites S into (L_strict + I + U, L_strict: strictly lower triangular, I: identity, U: upper triangular with nontrivial diagonal elements) which would save memory
    !this could be done by overwriting the non-sub-matrix part of S (one shell per iteration
    do cur_iter=1,n-1
       L(cur_iter,cur_iter)=1 !set diagonal element of L
       L(cur_iter+1:,cur_iter)=S(cur_iter+1:,1)/S(cur_iter,cur_iter)
       U(cur_iter,cur_iter:)=S(cur_iter,cur_iter:)
       S(cur_iter+1:,cur_iter+1:)=S(cur_iter+1:,cur_iter+1:)-rop(L(:,cur_iter),U(cur_iter,:)) !update Schur complement submatrix in S
    end do

    !may have to clean up last iteration manually
    !L(num_rows,num_rows)=1
    !U(num_rows,num_rows)=???
  end subroutine LU


end module decomp
