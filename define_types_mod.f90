module define_types

  use sizes
  
  implicit none

  type a_atmos
     logical::topdown
     real, dimension(nlayers):: temp,press,logP,dz
     double precision, dimension(nlayers)::opd_ext,opd_scat,gg
  end type a_atmos
  
  type a_gas
     logical::topdown
     integer:: id
     double precision, dimension(nlayers):: VMR
  end type a_gas


  type a_cloud
     logical::topdown
     integer::id
     double precision, dimension(nlayers)::density,rpeak,rsigma
  end type a_cloud

  save


  
end module define_types
