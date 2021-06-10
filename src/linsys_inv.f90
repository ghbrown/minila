
submodule (linsys) linsys_inv

contains
  module function inv(A) result(A_inv)
    !Computes the inverse of a square matrix
    implicit none
    real, intent(in) :: A(:,:)
    real, allocatable :: A_inv(:,:), id(:,:)
    integer :: j, n

    if (.not. is_square(A,n)) then
       print *, "ERROR: nonsquare matrices have no inverse"
       print *, "       (perhaps you'd like a psuedoinverse?)"
       return
    end if

    allocate(A_inv(n,n),id(n,n))

    !set identity matrix
    !move this to prim
    id=0.0
    do concurrent (j=1:n)
       id(j,j)=1.0
    end do

    A_inv=linsol(A,id) !solve A X = I using linsol

  end function inv

end submodule linsys_inv



  


