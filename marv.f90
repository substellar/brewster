subroutine marv(w1,w2,temp,logg,R2D2,ingasnum,logVMR,pcover,&
     do_clouds,cloudnum,cloudrad,cloudsig,cloudprof,&
     inlinetemps,inpress,inwavenum,inlinelist,cia,ciatemps,out_spec)

  use sizes
  use main
  

  !f2py integer, parameter :: nlayers
  !f2py integer, parameter :: nlinetemps
  !f2py integer, parameter :: ngas
  !f2py integer, parameter :: nclouds
  !f2py integer, parameter :: npatch
  !f2py integer, parameter :: nwave
  !f2py intent(in) logg,R2D2,pcover
  !f2py intent(in) w1,w2,temp,logVMR
  !f2py intent(in) ingasnum,do_clouds,cloudnum
  !f2py intent(in) cloudrad,cloudsig,cloudprof
  !f2py intent(in) inwavenum, inlinetemps,inpress
  !f2py intent(inout) cia, ciatemps
  !f2py intent(inout) inlinelist
  !f2py intent(out) out_spec

  real,dimension(4,nciatemps,nwave) :: cia
  real,dimension(nciatemps) :: ciatemps
  double precision,intent(inout) :: inlinelist(:,:,:,:)
  double precision,dimension(nlayers):: temp
  real :: R2D2,logg
  double precision :: w1,w2
  real,dimension(npatch)::pcover
  integer,dimension(npatch)::do_clouds
  character(len=10),dimension(ngas) :: gasname
  double precision,dimension(ngas) :: molmass
  integer,dimension(ngas)::ingasnum
  character(len=10),dimension(:),allocatable::gaslist,cloudlist
  double precision,dimension(:),allocatable:: masslist
  double precision,dimension(ngas,nlayers) :: logVMR
  integer,dimension(npatch,nclouds):: cloudnum
  character(len=10),dimension(npatch,nclouds) :: cloudname
  double precision,dimension(npatch,nlayers,nclouds) :: cloudrad
  double precision,dimension(npatch,nlayers,nclouds) :: cloudsig
  double precision,dimension(npatch,nlayers,nclouds) :: cloudprof
  double precision,dimension(2,nwave) :: out_spec
  double precision,dimension(nwave) :: inwavenum
  real,dimension(nlinetemps) :: inlinetemps
  real,dimension(nlayers) :: inpress
  integer:: maxgas,maxcloud,igas,icloud, idum1, idum2

  
  open(10,file="gaslist.dat", status='old')
  read(10,*) maxgas
  allocate(gaslist(maxgas), masslist(maxgas))
  do igas = 1, maxgas
     read(10,"(I4,A8,F9.5)") idum1,gaslist(igas),masslist(igas)
  end do
  close(10)


  
  do igas = 1, ngas
     gasname(igas) = adjustl(trim(gaslist(ingasnum(igas))))
     molmass(igas) = masslist(ingasnum(igas))
  end do

  
  open(10,file="cloudlist.dat", status='old')
  read(10,*) maxcloud
  allocate(cloudlist(maxcloud))
  do icloud = 1, maxcloud
     read(10,"(I3,A10)") idum2,cloudlist(icloud)
  end do
  close(10)

  
  do ipatch = 1, npatch
     do icloud = 1, nclouds
        cloudname(ipatch,icloud) = trim(adjustl(cloudlist(cloudnum(ipatch,icloud))))
     end do
  end do

  
  deallocate(cloudlist,gaslist,masslist)  


  call forward(w1,w2,temp,logg,R2D2,gasname,ingasnum,molmass,logVMR,pcover,&
       do_clouds,cloudname,cloudrad,cloudsig,cloudprof,&
       inlinetemps,inpress,inwavenum,inlinelist,cia,ciatemps,out_spec)

 
end subroutine marv


