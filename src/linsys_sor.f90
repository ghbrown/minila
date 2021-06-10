
submodule (linsys) linsys_sor

contains

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

end submodule linsys_sor

