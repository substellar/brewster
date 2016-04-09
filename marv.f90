subroutine marv(temp,logg,R2D2,ingasnum,logVMR,pcover,&
     do_clouds,incloudnum,cloudrad,cloudsig,cloudprof,&
     inlinetemps,inpress,inwavenum,inlinelist,cia,ciatemps,use_disort,outspec)

  use sizes
  use main
  

  !f2py integer, parameter :: nlayers
  !f2py integer, parameter :: nlinetemps
  !f2py integer, parameter :: ngas
  !f2py integer, parameter :: nclouds
  !f2py integer, parameter :: npatch
  !f2py intent(in) logg,R2D2
  !f2py intent(in) temp,logVMR
  !f2py intent(in) use_disort
  !f2py intent(in) inlinetemps,inpress
  !f2py intent(inout) cloudrad,cloudsig,cloudprof
  !f2py intent(inout) cia, ciatemps
  !f2py intent(inout) inlinelist, inwavenum
  !f2py intent(inout) do_clouds,pcover,ingasnum,incloudnum
  !f2py intent(out) out_spec

  real,intent(inout) :: cia(:,:,:)
  real,dimension(nciatemps) :: ciatemps
  double precision,intent(inout) :: inlinelist(:,:,:,:)
  double precision,dimension(nlayers):: temp
  real :: R2D2,logg
!  double precision :: w1,w2
  real,intent(inout) :: pcover(:)
  integer,intent(inout) ::do_clouds(:)
  character(len=10),dimension(:),allocatable :: gasname
  double precision,dimension(:),allocatable :: molmass
  integer,intent(inout) :: ingasnum(:)
  character(len=10),dimension(:),allocatable::gaslist,cloudlist
  double precision,dimension(:),allocatable:: masslist
  double precision,intent(inout) :: logVMR(:,:)
  integer,intent(inout) :: incloudnum(:,:)
  character(len=10),dimension(:,:),allocatable :: cloudname
  double precision,intent(inout) :: cloudrad(:,:,:)
  double precision,intent(inout) :: cloudsig(:,:,:)
  double precision,intent(inout) :: cloudprof(:,:,:)
  double precision,dimension(2,maxwave),intent(OUT):: outspec
  double precision,dimension(:,:),allocatable :: out_spec
  double precision,intent(inout) :: inwavenum(:)
  real,dimension(nlinetemps) :: inlinetemps
  real,dimension(nlayers) :: inpress
  integer:: maxgas,maxcloud,igas,icloud, idum1, idum2,use_disort

  call initwave(size(inwavenum))
  call initgas(size(ingasnum))
  call initpatch(size(do_clouds))

  allocate(molmass(ngas),gasname(ngas))
  allocate(cloudname(npatch,nclouds))
  
  
  allocate(out_spec(2,nwave))
  
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
        ! check if we're doing a specific cloud or a generic/mixed cloud
        if (incloudnum(ipatch,icloud) .eq. 99) then
           cloudname(ipatch,icloud) = "mixto"
        else
           cloudname(ipatch,icloud) = trim(adjustl(cloudlist(incloudnum(ipatch,icloud))))
        endif
     end do
  end do

  
  deallocate(cloudlist,gaslist,masslist)  
  call forward(temp,logg,R2D2,gasname,ingasnum,molmass,logVMR,pcover,&
       do_clouds,incloudnum,cloudname,cloudrad,cloudsig,cloudprof,&
       inlinetemps,inpress,inwavenum,inlinelist,cia,ciatemps,use_disort,out_spec)

  outspec = 0.0
  outspec(1,1:nwave)  = out_spec(1,1:nwave)
  outspec(2,1:nwave) = out_spec(2,1:nwave)

  

  deallocate(out_spec)
 
end subroutine marv


