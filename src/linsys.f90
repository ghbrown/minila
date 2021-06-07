
module linsys
  !Contains routines for the solution of linear systems, and computation of inverse.
  use prim
  use decomp
  implicit none

contains



  function linsol(A,b) result(x)
    !Solves the linear system Ax=b for x.
    !---Inputs---
    !A: full rank square matrix (these are checked)
    !b: right hand side of linear system
    !---Outputs---
    !x: solution to linear system
    real, intent(in) :: A(:,:)
    real, intent(in) :: b(size(A,2))
    real, allocatable :: x(:)
    real, allocatable, dimension(:,:) :: L, U
    real, allocatable :: b_permuted(:), y(:)
    integer, allocatable :: p(:) !row pivot vectors (such that PA=A(p,:)
    integer :: i, n

    if (.not. is_square(A,n)) then
       print *, "ERROR: cannot solve nonsquare linear system"
       return
    end if

    allocate(L(n,n), U(n,n), b_permuted(n), y(n), p(n))

    call lupp(A,L,U,p)

    if (minval(abs(get_diag(U))) == 0.0) then
       print *, "ERROR: linear system is singular"
       return
    end if
  
    b_permuted=b(p)
    y=fwd_sub(L,b_permuted) !solve L(Ux)=b_permuted for y=Ux
    x=bck_sub(U,y) !solve Ux=y=(Ux) for x

  end function linsol


  function inv(A) result(A_inv)
    !Computes the inverse of a square matrix
    implicit none
    real, intent(in) :: A(:,:)
    real, allocatable :: A_inv(:,:)
    integer :: n

    if (.not. is_square(A,n)) then
       print *, "ERROR: nonsquare matrices have no inverse"
       print *, "       (perhaps you'd like a psuedoinverse?)"
       return
    end if

    allocate(A_inv(n,n))

  end function inv

  !ITERATIVE METHODS
  function sor(A,b,rel_tol,max_iter) result(x)
    !computes solution to iterative system Ax=b via
    !successive over-relaxation (SOR)
    !WARNING: this implementation is (especially) cobbled together
    !---Inputs---
    !A: square matrix defining linear system coefficients
    !b: vector defining right hand side of linear system
    !rel_tol: relative tolerance of solution defining convergence
    !max_iter: maximum number of iterations defining convergence
    !---Outputs---
    !x: solution vector to Ax=b
    implicit none
    real, intent(in) :: A(:,:), rel_tol
    real, intent(in) :: b(:)
    integer, intent(in) :: max_iter
    real, allocatable :: x(:)
    real, allocatable :: x_prev(:), d(:), rhs(:)
    real, allocatable :: dpol(:,:), oupom1d(:,:)
    real :: omega, rel_change !taken from an of 3 cases in 
    integer :: i, n, cur_iter


    if (.not. is_square(A,n)) then
       print *, "ERROR: cannot solve nonsquare linear system"
       return
    end if

    allocate(x(n), x_prev(n), d(n), rhs(n), dpol(n,n), oupom1d(n,n))
    call random_number(x) !set random starting point for solution
    call random_number(x_prev) !set random previous vector


    omega=1.74 !relaxation parameter, taken as average of three cases
    ! in J. K. Reid's paper on subject of optimal relaxation parameter
    ! this does not mean it is a "good" value for any given problem,
    ! but it does mean it's reasonable
    d=get_diag(A) !diagonal of A (as vector)
    !form matrix  D + omega L
    dpol=set_diag(d)+omega*get_triang(A,"l")
    !form matrix  omega U + (omega - 1) D
    oupom1d= omega*get_triang(A,"u")+set_diag((omega-1)*d)

    rel_change=1.0e10 !fictitiously large number to enter while loop
    cur_iter=1

    do while ((rel_change > rel_tol) .and. (cur_iter <= max_iter))
       rhs=omega*b-matmul(oupom1d,x) !compute right hand side (vector) of over-relaxation equation
       x=fwd_sub(dpol,rhs) !update
       rel_change=norm2(x-x_prev)/norm2(x_prev) !relative change in solution
       x_prev=x
       cur_iter=cur_iter+1
    end do

  end function sor

  


end module linsys
