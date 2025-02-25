subroutine marv(temp,logg,R2D2,ingasname,molmass,logVMR,pcover,&
     do_clouds,incloudnum,cloudrad,cloudsig,cloudprof,&
     inlinetemps,inpress,inwavenum,inlinelist,cia,ciatemps,&
     use_disort,make_cl_pspec,make_oth_pspec,make_cf,do_bff,bff,outspec,&
     cl_phot_press,oth_phot_press,cfunc)

  use sizes
  use main
  

  !f2py integer, parameter :: nlinetemps
  !f2py intent(in) logg,R2D2,gasname
  !f2py intent(inout) temp,logVMR,inpress
  !f2py intent(in) use_disort,make_cl_pspec,make_oth_pspec,make_cf,do_bff
  !f2py intent(inout) inlinetemps
  !f2py intent(inout) cloudrad,cloudsig,cloudprof
  !f2py intent(inout) cia, ciatemps
  !f2py intent(inout) inlinelist, inwavenum,bff
  !f2py intent(inout) do_clouds,pcover,molmass,incloudnum
  !f2py intent(out) out_spec, cl_phot_press,oth_phot_press,cfunc

  real,intent(inout) :: cia(:,:,:)
  real,dimension(nciatemps) :: ciatemps
  double precision,intent(inout) :: inlinelist(:,:,:,:)
  double precision,intent(inout):: temp(:)
  real :: R2D2,logg
  double precision,intent(inout):: bff(:,:)
  real,intent(inout) :: pcover(:),molmass(:)
  integer,intent(inout) ::do_clouds(:)
  character(len=15),intent(in) :: ingasname(:)
  character(len=15),dimension(:),allocatable::cloudlist,gasname
  double precision,intent(inout) :: logVMR(:,:)
  integer,intent(inout) :: incloudnum(:,:)
  character(len=15),dimension(:,:),allocatable :: cloudname
  double precision,intent(inout) :: cloudrad(:,:,:)
  double precision,intent(inout) :: cloudsig(:,:,:)
  double precision,intent(inout) :: cloudprof(:,:,:)
  double precision,dimension(2,maxwave),intent(OUT):: outspec
  double precision,dimension(maxpatch,maxwave),intent(OUT):: cl_phot_press,oth_phot_press
  double precision,dimension(:,:),allocatable :: out_spec, clphotspec,othphotspec
  double precision,dimension(:,:,:),allocatable :: cf
  double precision,dimension(maxpatch,maxwave,maxlayers) :: cfunc
  double precision,intent(inout) :: inwavenum(:)
  !real,dimension(nlinetemps) :: inlinetemps
  real,intent(inout) :: inlinetemps(:)
  real,intent(inout) :: inpress(:)
  integer :: maxgas,maxcloud,igas,icloud, idum1, idum2
  integer :: use_disort,make_oth_pspec,make_cl_pspec,do_bff,make_cf
  logical :: othphot,clphot,do_cf
  
  call initlayers(size(inpress))
  call initwave(size(inwavenum))
  call initgas(size(molmass))
  call initpatch(size(do_clouds))
  call initcloud(size(cloudprof(1,1,:)))
  call inittemps(size(inlinetemps))

! allocate(molmass(ngas),gasname(ngas))
  allocate(cloudname(npatch,nclouds))

  allocate(out_spec(2,nwave))

  clphot = make_cl_pspec
  othphot = make_oth_pspec
  do_cf = make_cf
  
  allocate(gasname(size(ingasname)))
  gasname = ingasname
 
  open(10,file="data/cloudlist.dat", status='old')
  read(10,*) maxcloud
  allocate(cloudlist(maxcloud))
  do icloud = 1, maxcloud
     read(10,"(I3,A15)") idum2,cloudlist(icloud)
  end do
  close(10)
  
   do ipatch = 1, npatch
     do icloud = 1, nclouds
        ! check if we're doing a specific cloud or a generic/mixed cloud
        if (incloudnum(ipatch,icloud) .gt. 50) then
           cloudname(ipatch,icloud) = "mixto"
        else
           cloudname(ipatch,icloud) = trim(adjustl(cloudlist(incloudnum(ipatch,icloud))))
        endif
     end do
  end do

  
  deallocate(cloudlist)  
  call forward(temp,logg,R2D2,gasname,molmass,logVMR,pcover,&
       do_clouds,incloudnum,cloudname,cloudrad,cloudsig,cloudprof,&
       inlinetemps,inpress,inwavenum,inlinelist,cia,ciatemps,use_disort,&
       clphot,othphot,do_cf,do_bff,bff,out_spec,clphotspec,othphotspec,cf)

  outspec = 0.d0
  outspec(1,1:nwave)  = out_spec(1,1:nwave)
  outspec(2,1:nwave) = out_spec(2,1:nwave)

  cl_phot_press  = 0.d0
  if (clphot) then
     do ipatch = 1, npatch
        cl_phot_press(ipatch,1:nwave) = clphotspec(ipatch,1:nwave)
     end do
  end if
  oth_phot_press  = 0.d0
  if (othphot) then
     do ipatch = 1, npatch
        oth_phot_press(ipatch,1:nwave) = othphotspec(ipatch,1:nwave)
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

  deallocate(out_spec,clphotspec,othphotspec,cf)
 
end subroutine marv


