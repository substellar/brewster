 module setup_RT

contains

  subroutine run_RT(spectrum,clphotspec,othphotspec,cf,disorting,clphot,othphot,do_cf)

    use sizes
    use common_arrays
    use define_types
    use phys_const
    use atmos_ops

    
    implicit none

    integer, parameter :: MAXMOM = 16
    integer, parameter :: MAXPHI = 3
    integer, parameter :: MAXULV = 1
    integer, parameter :: MAXUMU = 16
    integer :: nlevel
    integer :: MAXCLY


    ! now set up where we want the fluxes.
    ! we just want the top of atmosphere to start
    ! can worry about bonus bits for tracking contributions later
    
    integer ::ntau
    double precision,dimension(MAXULV) :: utau
    double precision, allocatable, dimension(:,:) :: tau_cloud, tau_others   
    double precision,allocatable,dimension(:), INTENT(OUT):: spectrum
    double precision,allocatable,dimension(:,:), INTENT(OUT):: clphotspec,othphotspec
   
    double precision,allocatable,dimension(:,:,:), INTENT(OUT):: cf
    double precision, dimension(nlayers) :: DTAUC, SSALB, COSBAR
    double precision, dimension(nlayers+1) :: temper
    double precision :: WVNMLO, WVNMHI, wint, tau1, tau2, p1, p2
    double precision :: taup_cl,taup_oth,tau
    integer :: ipatch,ilayer, iwave
    double precision,dimension(0:MAXMOM,nlayers) :: PMOM
    double precision:: phi(maxphi),phi0, umu(maxumu), umu0
    double precision,dimension(MAXULV) :: RFLDIR,RFLDN,FLUP,DFDT,UAVG
    double precision,dimension(nlayers+1) :: gflup,fdi
    double precision,dimension(maxumu,maxulv,maxphi) :: UU
    double precision,dimension(maxumu) :: ALBMED, TRNMED
    double precision:: BTEMP, TTEMP
    double precision,allocatable,dimension(:):: upflux
    logical, dimension(5) :: PRNT
    character(len=0):: HEADER
    integer :: NUMU,NSTR,NMOM,NLYR, NPHI, IBCND
    logical :: LAMBER,  PLANK,  USRTAU, USRANG, ONLYFL,disorting
    logical :: clphot,othphot,othdone, cldone,do_cf
    double precision :: FBEAM, FISOT,  ALBEDO , ACCUR, TEMIS
    double precision :: bbplk
    external :: bbplk

    nlevel = nlayers+1
    MAXCLY = nlayers
!    integer :: nw1, nw2
    
    HEADER = ""

    PRNT = [.FALSE., .FALSE.,.FALSE.,.FALSE.,.FALSE.]

    ! set values for non-parameter disort input
    NSTR = 8
    NMOM = 8
    NUMU = 8 ! same as NSTR
    NLYR = nlayers
    NPHI = 0
    NTAU = 1
    UTAU = 0.0
    
    USRTAU    = .TRUE.  ! return quantities at each pre-defined layer
    USRANG    = .FALSE.
    ONLYFL    = .TRUE.   ! return only fluxes      
    IBCND     = 0        ! boundary conditions, change if Planck is false
    FBEAM     = 0.0      ! parallel beam at the top. 0 for stars
    ! if this is not 0, need to specify umu0 and phi0
    umu0 = 0.6667
    FISOT     = 0.0      ! not sure, was 1.0 / PI
    LAMBER    = .TRUE.   ! don't care about lower boundary
    ALBEDO = 0.000000001
    ACCUR = 0.01 
    PLANK    = .TRUE.    ! need this to use temperature structure
    TEMIS = 0.0        ! need to give top layer a bit of emissivity 
    

    ! all temps the same across patches
    call set_temp_levels(temper)

    allocate(upflux(nwave),spectrum(nwave))
    allocate(clphotspec(npatch,nwave),othphotspec(npatch,nwave))
    allocate(cf(npatch,nwave,nlayers))
    allocate(tau_cloud(nlayers,nwave),tau_others(nlayers,nwave))

    clphotspec = 0d0
    othphotspec = 0d0
    spectrum = 0d0
    
    do ipatch = 1, npatch
       ! add up the taus to get extinction
       do ilayer = 1, nlayers

          tau_cloud(ilayer,:) = patch(ipatch)%atm(ilayer)%opd_ext

          tau_others(ilayer,:) = patch(ipatch)%atm(ilayer)%opd_lines + &
               patch(ipatch)%atm(ilayer)%opd_CIA + &
               patch(ipatch)%atm(ilayer)%opd_rayl + &
               patch(ipatch)%atm(ilayer)%opd_hmbff
          
          
          patch(ipatch)%atm(ilayer)%opd_ext = &
               patch(ipatch)%atm(ilayer)%opd_ext + &
               patch(ipatch)%atm(ilayer)%opd_lines + &
               patch(ipatch)%atm(ilayer)%opd_CIA + &
               patch(ipatch)%atm(ilayer)%opd_rayl + &
               patch(ipatch)%atm(ilayer)%opd_hmbff

          patch(ipatch)%atm(ilayer)%opd_scat = &
               patch(ipatch)%atm(ilayer)%opd_scat + &
               patch(ipatch)%atm(ilayer)%opd_rayl

                    
       end do
     



       upflux = 0.0
       !$OMP PARALLEL DO num_threads(8) default(SHARED) PRIVATE(iwave,ilayer,cldone,othdone,tau1,tau2,p1,p2,WVNMLO,DTAUC,SSALB,COSBAR,ALBEDO,gflup,fdi)
       do iwave = 1, nwave
       
       ! need PMOM
        ! get this calling Mark's GETMOM code
        ! 6 is Rayleigh+HG
        ! 2 ia just Rayleigh
          do ilayer = 1, nlayers

             if (patch(ipatch)%cloudy .ne. 0) then
                SSALB(ilayer) = patch(ipatch)%atm(ilayer)%opd_scat(iwave) / &
                     patch(ipatch)%atm(ilayer)%opd_ext(iwave)             
                ! test lines
                !SSALB(ilayer) = 0.5d0
                !patch(ipatch)%atm(ilayer)%gg(iwave) = 0.0d0
                !if (disorting) call &
                     ! GETMOM(6,patch(ipatch)%atm(ilayer)%gg(iwave),&
                     ! NMOM,0.0, PMOM(0,ilayer))
                COSBAR(ilayer) = patch(ipatch)%atm(ilayer)%gg(iwave)
             else     
                ! CALL GETMOM( 2, 0.0, NMOM, 0.0, PMOM(0:nmom,ilayer))
                SSALB(ilayer) = patch(ipatch)%atm(ilayer)%opd_scat(iwave) / &
                     patch(ipatch)%atm(ilayer)%opd_ext(iwave)             
                COSBAR(ilayer) = 0.0
             end if
             
             ! TEST LINE
             !DTAUC(ilayer) =  patch(ipatch)%atm(ilayer)%dp /1000.

          enddo ! layer loop

!          write(*,*) sum(dtauc)
          DTAUC = 0.d0

          ! get the pressure level where tau_cloud = taup_cl
          ! set reference tau for cloud taup_cl
          taup_cl = 1.0
          ! set reference tau for others taup_oth
          taup_oth = 1.0


          cldone = .false.
          othdone = .false.

          ! Run through the layers...
          
          do ilayer = 1, nlayers
             
             !put optical depth into the right variable for radtran
             DTAUC(ilayer) = patch(ipatch)%atm(ilayer)%opd_ext(iwave)

             ! now sort out the diagnostics for the photospheres
             if (clphot .and. .not. (cldone)) then
                
                ! this bit calculates the pressure level where tau_cloud
                ! reaches some value set above. Activate in python code with
                ! gnostics
                if (sum(tau_cloud(1:ilayer,iwave)) .gt. taup_cl) then
                   cldone = .true.
                   tau2 = sum(tau_cloud(1:ilayer,iwave))
                   tau1 = tau2 - tau_cloud(ilayer,iwave)
                   
                   
                   if (ilayer .eq. nlayers) then
                      p1 = exp((0.5)*(log(patch(ipatch)%atm(ilayer-1)%press * patch(ipatch)%atm(ilayer)%press)))
                   else
                      p1 = exp(((1.5)*log(patch(ipatch)%atm(ilayer)%press)) - &
                           ((0.5)*log(patch(ipatch)%atm(ilayer+1)%press)))
                   end if
                   
                   clphotspec(ipatch,iwave) = p1 +((taup_cl - tau1) * &
                        patch(ipatch)%atm(ilayer)%dp / tau_cloud(ilayer,iwave))
                end if
             end if
             
             if (othphot .and. .not. (othdone)) then
                ! this bit calculates the pressure level where tau_other
                ! (i.e. not clouds) reaches some value set above.
                ! Activate in python code with gnostics
                
                if (sum(tau_others(1:ilayer,iwave)) .gt. taup_oth) then
                   othdone = .true.
                   tau2 = sum(tau_others(1:ilayer,iwave))
                   tau1 = tau2 - tau_others(ilayer,iwave)
                   
                   
                   if (ilayer .eq. nlayers) then
                      p1 = exp((0.5)*(log(patch(ipatch)%atm(ilayer-1)%press * patch(ipatch)%atm(ilayer)%press)))
                   else
                      p1 = exp(((1.5)*log(patch(ipatch)%atm(ilayer)%press)) - &
                           ((0.5)*log(patch(ipatch)%atm(ilayer+1)%press)))
                   end if
                   
                   othphotspec(ipatch,iwave) = p1 +((taup_oth - tau1) * &
                        patch(ipatch)%atm(ilayer)%dp / tau_others(ilayer,iwave))
                end if
             end if
             
             
          end do
          
          
          if (disorting) then 
             
             
             ! set up wavenumber interval....
             
             
             if (iwave .eq. 1) then 
                WVNMLO = wavenum(1) - 0.5*(wavenum(2) - wavenum(1))
                WVNMHI =  wavenum(1) + 0.5*(wavenum(2) - wavenum(1))
             else if (iwave .eq. nwave) then
                WVNMLO = wavenum(nwave) - 0.5*(wavenum(nwave) - wavenum(nwave-1))
                WVNMHI =  wavenum(nwave) + 0.5*(wavenum(nwave) - wavenum(nwave-1))
             else
                WVNMLO = wavenum(iwave) - 0.5*(wavenum(iwave) - wavenum(iwave-1))
                WVNMHI = wavenum(iwave) + 0.5*(wavenum(iwave+1) -wavenum(iwave))
             end if
             
             ! set top and bottom boundary temperatures
             BTEMP = temper(nlayers+1)
             TTEMP = temper(1)
             
             ! convert to flux density W/m2/um
             ! need interval in um not cm^-1
             wint = (1.0e4 / wavenum(iwave)) * ((WVNMHI - WVNMLO)/ wavenum(iwave))
             
             call DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO, &
                  WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, &
                  UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, &
                  FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS,&
                  PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY,&
                  MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, RFLDN,&
                  FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )
             
             upflux(iwave) = FLUP(1) / wint
             
          else
             call gfluxi(temper,DTAUC,SSALB,COSBAR,wavenum(iwave),ALBEDO,gflup,&
                  fdi)
             upflux(iwave) = gflup(1) !/ wint
          endif
          
       end do ! wave loop
       !$OMP END PARALLEL DO
       spectrum = spectrum + (upflux*patch(ipatch)%cover)
       
    end do ! patch loop
    
       
    deallocate(upflux)
    ! Get the contribution function for them that want it...
    if (do_CF) then
       do ipatch = 1, npatch
          do iwave = 1, nwave
             tau = 0.d0
             do ilayer = 1, nlayers
                tau = tau + patch(ipatch)%atm(ilayer)%opd_ext(iwave)
                cf(ipatch,iwave,ilayer)  = &
                     bbplk(wavenum(iwave),patch(ipatch)%atm(ilayer)%temp) * &
                     patch(ipatch)%atm(ilayer)%opd_ext(iwave) &
                     / exp(tau)
             end do
          end do
       end do
    end if

    
  end subroutine run_RT

    
  
end module setup_RT
