module decomp
  !Contains subroutines for matrix decompositions: LU 
  !and eventually: QR, eigenvalue
  implicit none

  use primitives

contains

  subroutine LU(A,L,U)
    implicit none
    real, intent(in) :: A(:,:)
    integer :: num_rows=size(A,1)
    integer :: num_cols=size(A,2)
    real, intent(out) :: L(num_rows,num_rows), U(num_rows,num_rows)
    real :: S(num_rows,num_rows)
    if (num_rows /= num_cols) then !check if A is square
       print *, 'ERROR: LU is invalid for nonsquare matrices'
       return !LU is not valid for nonsquare matrices, exit early
    end if

    S=A
    L=0.0
    U=0.0
    !in principle it should be possible to do this whole algorithm without recursion by maintining an "A-sized" S throughout, but only using progressively less and less of it (shrinking the used submatrix) by also indexing it with cur_iter
    !this will take more memory and about twice as long, but it's infinitely more readable
    do cur_iter=1,n-1
       L(cur_iter,cur_iter)=1 !set diagonal element of L
       L(cur_iter+1:,cur_iter)=S(cur_iter+1:,1)/S(cur_iter,cur_iter)
       U(cur_iter,cur_iter:)=S(cur_iter,cur_iter:)
       S=S-r1p(L(:,cur_iter),U(cur_iter,:)) 
    end do
    !may have to clean up last iteration manually
    !L(num_rows,num_rows)=1
    !U(num_rows,num_rows)=???
  end subroutine LU


end module decomp
