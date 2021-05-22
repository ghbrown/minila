
program test_primitives
  use primitives
  implicit none

  real, dimension(3) :: rop_vec1=(/1,2,3/)
  real, dimension(3) :: rop_vec2=(/4,5,6/)
  real, dimension(3,3) :: rop_mat_ans=transpose(reshape((/4.,5.,6.,8.,10.,12.,12.,15.,18./),(/3,3/)))
  real :: rop_res_mat(3,3) !rank one product residual matrix

  rop_res_mat=abs(rop(rop_vec1,rop_vec2)-rop_mat_ans)
  if ((sum(rop_res_mat)/size(rop_res_mat)) > 1e-10) then
     print *, "Rank one product (rop) function failed."
  end if

end program test_primitives

