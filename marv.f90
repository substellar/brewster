subroutine marv(w1,w2,temp,logg,R2D2,gasnum,VMR,pcover,&
     do_clouds,cloudnum,cloudrad,cloudsig,cloudprof,out_spec)

  use sizes
  use main

  !f2py integer, parameter :: nlayers
  !f2py integer, parameter :: ngas
  !f2py integer, parameter :: nclouds
  !f2py integer, parameter :: npatch
  !f2py integer, parameter :: nwave
   
  real,dimension(nlayers)::temp
  real ::R2D2,logg,w1,w2
  real,dimension(npatch) :: pcover
  integer,dimension(npatch) :: do_clouds
  integer,dimension(ngas)::gasnum
  character(len=10),dimension(ngas) :: gasname
  double precision, dimension(ngas) :: molmass
  character(len=10),dimension(:),allocatable::gaslist,cloudlist
  double precision, dimension(ngas,nlayers) :: VMR
  integer,dimension(npatch,nclouds):: cloudnum
  character(len=10),dimension(npatch,nclouds) :: cloudname
  double precision, dimension(npatch,nlayers,nclouds) :: cloudrad
  double precision, dimension(npatch,nlayers,nclouds) :: cloudsig
  double precision, dimension(npatch,nlayers,nclouds) :: cloudprof
  real,dimension(2,nwave) :: out_spec
  integer:: maxgas,maxcloud,igas,icloud, idum1, idum2
  !f2py intent(in) w1,w2,temp,logg,R2D2,gasnum,VMR,pcover,do_clouds
  !f2py intent(in) cloudnum,cloudrad,cloudsig,cloudprof
  !f2py intent(out) out_spec
  !f2py depend(nlayers,npatch,nclouds) cloudrad,cloudsig,cloudprof
  !f2py depend(ngas,nlayers) VMR
  !f2py depend(npatch,nclouds) cloudnum,cloudname
  !f2py depend(npatch) pcover,do_clouds 
  !f2py depend(ngas) gasnum
  !f2py depend(nlayers) temp
  !f2py depend(nwave) out_spec

  
  write(*,*) temp
  
  open(10,file="gaslist.dat", status='old')
  read(10,*) maxgas
  allocate(gaslist(maxgas))
  do igas = 1, maxgas
     read(10,"(I4,A8,F9.5)") idum1,gaslist(igas),molmass(igas)
  end do
  close(10)
  
  do igas = 1, ngas
     gasname(igas) = adjustl(trim(gaslist(gasnum(igas))))
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
  
  deallocate(cloudlist,gaslist)  
  
  call forward(w1,w2,temp,logg,R2D2,gasname,molmass,VMR,pcover,&
       do_clouds,cloudname,cloudrad,cloudsig,cloudprof,out_spec)



end subroutine marv
