module decomp
implicit none

contains

  subroutine LU(A,L,U)
    implicit none
    real, intent(in) :: A(:,:)
    integer :: num_rows=size(A,1)
    integer :: num_cols=size(A,2)
    if (num_rows /= num_cols) then !check if A is square
       print *, 'ERROR: LU is invalid for nonsquare matrices'
       return !LU is not valid for nonsquare matrices, exit early
    end if
    real, intent(out) :: L(num_rows,num_rows), U(num_rows,num_rows)
    real :: S(num_rows,num_rows)

    S=A
    L=0.0
    U=0.0
    !this principle is that it should be possible to do this whole algorithm without recursion by maintining an "A-sized" S throughout, but only using progressively less and less of it (shrinking the submatrix) by also indexing it with cur_iter
    !this will take more memory and about twice as long, but it's infinitely more readable
    do cur_iter=1,n-1
       L(cur_iter,cur_iter)=1 !set diagonal element of L
       L(cur_iter+1:,cur_iter)=S(cur_iter+1:,1)/S(cur_iter,cur_iter)
       U(cur_iter,cur_iter:)=S(cur_iter,cur_iter:)
       S=S-matmul(L(:,cur_iter),U(cur_iter,:) 
    end do
    !may have to clean up last iteration manually
    !L(num_rows,num_rows)=1
    !U(num_rows,num_rows)=???
  end subroutine

  subroutine QR(A,Q,R)
    implicit none
    real, intent(in) :: A
    real, intent(out) :: Q, R
  end subroutine

end module decomp
