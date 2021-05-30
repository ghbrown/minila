
module linsys
  !linear system solvers
  !for square linear system, check for singular A by looking at the diagonla elements of matrix U
  !if A is singular, one of U's diagonal elements will be 0

  !check if A is square
  !call lupp(A,L,U,P)
  !if (min(abs(get_diag(U))) < sing_tol) then !check diagonal of u to see if A is invertible
  !print *, "ERROR: Attempted to solve singular linear system"
  !else
  !b=b(P) !pivot swap b using P vector
  !y=fwdsub(L,b)
  !x=bcksub(U,y)

end module linsys
