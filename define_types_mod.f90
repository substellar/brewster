module define_types

  use sizes
  
  implicit none

  type a_gas
     character(len=10):: name
     integer :: num
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
     double precision:: temp
     double precision :: press,logP,dz,dp,ndens,fH2,fHe,mu
     double precision, allocatable,dimension(:) ::opd_ext,opd_scat,gg
     double precision, allocatable,dimension(:) :: opd_lines,opd_CIA,opd_rayl
     type(a_gas),allocatable,dimension(:) :: gas
     type(a_cloud) :: cloud(nclouds)
  end type a_layer

  type a_patch
     integer:: index
     logical :: cloudy
     real:: cover
     type(a_layer)::atm(nlayers)
  end type a_patch

  
  save

contains
  subroutine init_column(col)
    type(a_layer) :: col(nlayers)
    integer:: ipatch,ilayer
    do ilayer = 1, nlayers

       if ( .NOT. allocated (col(ilayer)%gas)) &
            allocate (col(ilayer)%gas(ngas))
       if ( .NOT. allocated (col(ilayer)%opd_ext)) &
            allocate (col(ilayer)%opd_ext(nwave))
       if ( .NOT. allocated (col(ilayer)%opd_scat)) &
            allocate (col(ilayer)%opd_scat(nwave))
       if ( .NOT. allocated (col(ilayer)%opd_lines)) &
            allocate (col(ilayer)%opd_lines(nwave))
       if ( .NOT. allocated (col(ilayer)%opd_CIA)) &
            allocate (col(ilayer)%opd_CIA(nwave))
       if ( .NOT. allocated (col(ilayer)%opd_rayl)) &
            allocate (col(ilayer)%opd_rayl(nwave))
       if ( .NOT. allocated (col(ilayer)%gg)) &
            allocate (col(ilayer)%gg(nwave))

       
    end do
  end subroutine init_column

end module define_types
