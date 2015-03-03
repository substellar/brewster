module common_arrays_mod

  implicit none

  ! declares arrays that will be used alot

  integer, dimension(ngas) :: gasid
  real, dimension(nlayers):: temp,dz,press,logP
  double precision, dimension(ngas,nlayers):: VMR 


  ! we'll also set the values of the pressure array here



  
  press=[1.e-02, 3.e-02, 1.e-1, 3.e-1,1.e0, 3.e0,1.e1, 3.e1, 1.e2, &
           3.e2,1.e3, 3.e3,1.e4, 3.e4,1.e5,3.e5] 

  logP = log10(press)




  save


end module common_arrays_mod
