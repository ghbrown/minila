
submodule (linsys) linsys_inverse

contains

  module function inverse(A) result(A_inv)
    !Computes the inverse of a square matrix
    implicit none
    real, intent(in) :: A(:,:)
    real, allocatable :: A_inv(:,:), Id(:,:)
    integer :: j, n

    if (.not. is_square(A,n)) then
       print *, "ERROR: nonsquare matrices have no inverse"
       print *, "       (perhaps you'd like a psuedoinverse?)"
       return
    end if

    allocate(A_inv(n,n),Id(n,n))

    Id=identity(n) !make identity matrix of dimension n

    A_inv=linsol(A,Id) !solve A X = Id using linsol, giving X = A^(-1)
  end function inverse

end submodule linsys_inverse



  


