module define_types

  use sizes
  
  implicit none

    type a_gas
     character(len=10):: name
     double precision:: VMR
     double precision :: molmass
  end type a_gas


  type a_cloud
     character(len=10):: name
     double precision :: density,rg,rsig
  end type a_cloud


  type a_layer
     ! set an index to identify the layer in the fixed pressure scale
     ! layer 1 is top of atmosphere!!
     integer:: index
     real:: temp
     double precision :: press,logP,dz,ndens,fH2,fHe,mu
     double precision, dimension(nwave) ::opd_ext,opd_scat,gg,opd_lines,opd_CIA
     type(a_gas) :: gas(ngas)
     type(a_cloud) :: cloud(nclouds)
  end type a_layer
  

  type a_patch
     integer:: index
     logical :: cloudy
     real:: cover
     type(a_layer)::atm(nlayers)
  end type a_patch
  
  save


  
end module define_types
