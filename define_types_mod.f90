module define_types

  use sizes
  
  implicit none

    type a_gas
     character(len=10):: name
     double precision:: VMR
  end type a_gas


  type a_cloud
     character(len=10):: name
     double precision :: density,rpeak,rsigma
  end type a_cloud


  type a_layer
     ! set an index to identify the layer in the fixed pressure scale
     ! layer 1 is top of atmosphere!!
     integer:: index
     real:: temp
     double precision :: press,logP,dz,ndens
     double precision, dimension(nwave) ::opd_ext,opd_scat,gg,opd_lines,opd_CIA
     type(a_gas) :: gas(ngas)
     type(a_cloud) :: cloud(ncloud)
  end type a_layer
  

  save


  
end module define_types
