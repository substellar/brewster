subroutine marv(temp,logg,R2D2,ingasnum,logVMR,pcover,&
     do_clouds,incloudnum,cloudrad,cloudsig,cloudprof,&
     inlinetemps,inpress,inwavenum,inlinelist,cia,ciatemps,&
     use_disort,make_pspec,make_tspec,make_cf,do_bff,bff,outspec,&
     phot_press,tau_spec,cfunc)

  use sizes
  use main
  

  !f2py integer, parameter :: nlinetemps
  !f2py intent(in) logg,R2D2
  !f2py intent(inout) temp,logVMR,inpress
  !f2py intent(in) use_disort,make_pspec,make_tspec,make_cf,do_bff
  !f2py intent(in) inlinetemps
  !f2py intent(inout) cloudrad,cloudsig,cloudprof
  !f2py intent(inout) cia, ciatemps
  !f2py intent(inout) inlinelist, inwavenum,bff
  !f2py intent(inout) do_clouds,pcover,ingasnum,incloudnum
  !f2py intent(out) out_spec, phot_press,cfunc

  real,intent(inout) :: cia(:,:,:)
  real,dimension(nciatemps) :: ciatemps
  double precision,intent(inout) :: inlinelist(:,:,:,:)
  double precision,intent(inout):: temp(:)
  real :: R2D2,logg
  double precision,intent(inout):: bff(:,:)
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
  double precision,dimension(maxpatch,maxwave),intent(OUT):: phot_press,tau_spec
  double precision,dimension(:,:),allocatable :: out_spec, photspec,tauspec
  double precision,dimension(:,:,:),allocatable :: cf
  double precision,dimension(maxpatch,maxwave,maxlayers) :: cfunc
  double precision,intent(inout) :: inwavenum(:)
  real,dimension(nlinetemps) :: inlinetemps
  real,intent(inout) :: inpress(:)
  integer :: maxgas,maxcloud,igas,icloud, idum1, idum2
  integer :: use_disort,make_pspec,make_tspec,do_bff,make_cf
  logical :: pspec,tspec,do_cf
  
  call initlayers(size(inpress))
  call initwave(size(inwavenum))
  call initgas(size(ingasnum))
  call initpatch(size(do_clouds))
  call initcloud(size(cloudprof(1,1,:)))

  allocate(molmass(ngas),gasname(ngas))
  allocate(cloudname(npatch,nclouds))

  allocate(out_spec(2,nwave))

  pspec = make_pspec
  tspec = make_tspec
  do_cf = make_cf
  
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
       inlinetemps,inpress,inwavenum,inlinelist,cia,ciatemps,use_disort,pspec,tspec,do_cf,do_bff,bff,out_spec,photspec,tauspec,cf)

  outspec = 0.d0
  outspec(1,1:nwave)  = out_spec(1,1:nwave)
  outspec(2,1:nwave) = out_spec(2,1:nwave)

  phot_press  = 0.d0
  if (pspec) then
     do ipatch = 1, npatch
        phot_press(ipatch,1:nwave) = photspec(ipatch,1:nwave)
     end do
  end if
  tau_spec  = 0.d0
  if (tspec) then
     do ipatch = 1, npatch
        tau_spec(ipatch,1:nwave) = tauspec(ipatch,1:nwave)
     end do
  end if

  cfunc = 0.d0
  if (do_cf) then
     do ipatch = 1, npatch
        do ilayer = 1, nlayers
           cfunc(ipatch,1:nwave,ilayer)  = cf(ipatch,1:nwave,ilayer)
        end do
     end do
  end if

  deallocate(out_spec,photspec,tauspec,cf)
 
end subroutine marv


