
program test_prim
  implicit none

  real :: abs_tol=1e-10

  print *
  print *, "---PRIMITIVE TESTS---"
  call test_rop(abs_tol)
  call test_strang(abs_tol)
  print *

end program test_prim

subroutine test_rop(abs_tol)
  use prim

  implicit none
  real, intent(in) :: abs_tol
  real, dimension(3) :: rop_vec1=(/1,2,3/)
  real, dimension(3) :: rop_vec2=(/4,5,6/)
  real, dimension(3,3) :: rop_mat_ans=transpose(reshape((/4.,5.,6.,8.,10.,12.,12.,15.,18./),(/3,3/)))
  real :: rop_res_mat(3,3) !rank one product residual matrix

  rop_res_mat=abs(rop(rop_vec1,rop_vec2)-rop_mat_ans)
  if ((sum(rop_res_mat)/size(rop_res_mat)) > abs_tol) then
     print *, "  FAILED,    rop (rank one product)"
  else
     print *, "  PASSED,    rop (rank one product)"
  end if
end subroutine test_rop

subroutine test_strang(abs_tol)
  use prim

  implicit none
  real, intent(in) :: abs_tol
  real, dimension(3,3) :: T_upper=transpose(reshape((/1.,2.,3.,0.,4.,5.,0.,0.,6./),(/3,3/)))
  real, dimension(3) :: b=(/14.,23.,18./)
  real, dimension(3) :: x_ans=(/1.,2.,3./)
  real, dimension(3) :: x

  !test back substitution
  x=strang(T_upper,b)
  if (norm2(x-x_ans) > abs_tol) then
     print *, "  FAILED,    strang (triangular solve)" 
  else
     print *, "  PASSED,    strang (triangular solve)"
  end if

  !test forward substitution
  !x=strang(T_lower,b)
  

end subroutine test_strang






