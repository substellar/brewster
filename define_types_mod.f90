module define_types

  use sizes
  
  implicit none

  type a_gas
     character(len=10):: name
     double precision:: VMR
     double precision :: molmass
  end type a_gas


  type a_cloud
     character(len=15):: name
     double precision :: dtau1,rg,rsig
  end type a_cloud


  type a_layer
     ! set an index to identify the layer in the fixed pressure scale
     ! layer 1 is top of atmosphere!!
     integer:: index
     double precision:: temp
     double precision :: press,logP,dz,dp,ndens,fH2,fHe,fHmin,fH,fe,mu
     double precision, allocatable,dimension(:) ::opd_ext,opd_scat,gg
     double precision, allocatable,dimension(:) :: opd_lines,opd_CIA,opd_rayl,opd_hmbff
     type(a_gas),allocatable,dimension(:) :: gas
     type(a_cloud),allocatable,dimension(:) :: cloud
  end type a_layer

  type a_patch
     integer:: index
     integer :: cloudy
     real:: cover
     type(a_layer),allocatable,dimension(:):: atm
  end type a_patch

  
  save

contains
  subroutine init_column(col)
    type(a_layer),allocatable,dimension(:) :: col
    integer:: ipatch,ilayer
    if ( .NOT. allocated (col)) &
         allocate (col(nlayers))
    do ilayer = 1, nlayers
       if ( .NOT. allocated (col(ilayer)%gas)) &
            allocate (col(ilayer)%gas(ngas))
       if ( .NOT. allocated (col(ilayer)%cloud)) &
            allocate (col(ilayer)%cloud(nclouds))
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
       if ( .NOT. allocated (col(ilayer)%opd_hmbff)) &
            allocate (col(ilayer)%opd_hmbff(nwave))

       
    end do
  end subroutine init_column

end module define_types
