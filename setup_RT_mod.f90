 module setup_RT

contains

  subroutine run_RT(spectrum,photspec,tauspec,disorting,pspec,tspec)

    use sizes
    use common_arrays
    use define_types
    use phys_const
    use atmos_ops

    
    implicit none

    integer, parameter :: nlevel = nlayers+1
    integer, parameter :: MAXCLY = nlayers
    integer, parameter :: MAXMOM = 16
    integer, parameter :: MAXPHI = 3
    integer, parameter :: MAXULV = 1
    integer, parameter :: MAXUMU = 16


    ! now set up where we want the fluxes.
    ! we just want the top of atmosphere to start
    ! can worry about bonus bits for tracking contributions later
    
    integer ::ntau
    double precision,dimension(MAXULV) :: utau
    
    double precision,allocatable,dimension(:), INTENT(OUT):: spectrum
    double precision,allocatable,dimension(:,:), INTENT(OUT):: photspec,tauspec    
    double precision, dimension(nlayers) :: DTAUC, SSALB, COSBAR
    double precision, dimension(nlayers+1) :: temper
    double precision :: WVNMLO, WVNMHI, wint, tau1, tau2, p1, p2,ptau,taup
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
    logical :: pspec,tspec,tdone, pdone
    double precision :: FBEAM, FISOT,  ALBEDO , ACCUR, TEMIS
    


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
    allocate(photspec(npatch,nwave),tauspec(npatch,nwave))
    
    spectrum = 0.0
    
    do ipatch = 1, npatch
       ! add up the taus to get extinction
       do ilayer = 1, nlayers
        
          patch(ipatch)%atm(ilayer)%opd_ext = &
               patch(ipatch)%atm(ilayer)%opd_ext + &
               patch(ipatch)%atm(ilayer)%opd_lines + &
               patch(ipatch)%atm(ilayer)%opd_CIA + &
               patch(ipatch)%atm(ilayer)%opd_rayl

          patch(ipatch)%atm(ilayer)%opd_scat = &
               patch(ipatch)%atm(ilayer)%opd_scat + &
               patch(ipatch)%atm(ilayer)%opd_rayl
       end do
     



       upflux = 0.0
       do iwave = 1, nwave
       
       ! need PMOM
        ! get this calling Mark's GETMOM code
        ! 6 is Rayleigh+HG
        ! 2 ia just Rayleigh
          do ilayer = 1, nlayers

             if (patch(ipatch)%cloudy) then
                SSALB(ilayer) = patch(ipatch)%atm(ilayer)%opd_scat(iwave) / &
                     patch(ipatch)%atm(ilayer)%opd_ext(iwave)             
                ! test lines
                !SSALB(ilayer) = 0.5d0
                !patch(ipatch)%atm(ilayer)%gg(iwave) = 0.0d0
                if (disorting) call &
                     GETMOM(6,patch(ipatch)%atm(ilayer)%gg(iwave),&
                     NMOM,0.0, PMOM(0,ilayer))
                COSBAR(ilayer) = patch(ipatch)%atm(ilayer)%gg(iwave)
             else     
                CALL GETMOM( 2, 0.0, NMOM, 0.0, PMOM(0:nmom,ilayer))
                SSALB(ilayer) = patch(ipatch)%atm(ilayer)%opd_scat(iwave) / &
                     patch(ipatch)%atm(ilayer)%opd_ext(iwave)             
                COSBAR(ilayer) = 0.0
             end if
             
             ! TEST LINE
             !DTAUC(ilayer) =  patch(ipatch)%atm(ilayer)%dp /1000.

          enddo ! layer loop

!          write(*,*) sum(dtauc)
          DTAUC = 0.d0

          ! get the pressure level where tau = taup
          ! set taup
          taup = 3.0
          ! also  get tau at set pressure level
          ptau = 1.0

          pdone = .false.
          tdone = .false.
          do ilayer = 1, nlayers
             
             DTAUC(ilayer) = patch(ipatch)%atm(ilayer)%opd_ext(iwave)

             if (pspec .and. .not. (pdone)) then
                if (sum(dtauc) .gt. taup) then
                   pdone = .true.
                   tau2 = sum(dtauc)
                   tau1 = tau2 - dtauc(ilayer)

                   if (ilayer .eq. nlayers) then
                      p1 = exp((0.5)*(log(patch(ipatch)%atm(ilayer-1)%press * patch(ipatch)%atm(ilayer)%press)))
                   else
                      p1 = exp(((1.5)*log(patch(ipatch)%atm(ilayer)%press)) - &
               ((0.5)*log(patch(ipatch)%atm(ilayer+1)%press)))
                   end if

                   photspec(ipatch,iwave) = p1 +((taup - tau1) * &
                        patch(ipatch)%atm(ilayer)%dp / dtauc(ilayer))
                   
                end if
             end if

             if (tspec .and. .not. (tdone)) then
                if (sum(patch(ipatch)%atm(1:ilayer)%dp) .gt. ptau) then
                   tdone = .true.
                   p2 = sum(patch(ipatch)%atm(1:ilayer)%dp)
                   p1 = p2 - patch(ipatch)%atm(ilayer)%dp

                   tau1 = sum(dtauc) - dtauc(ilayer)

                   tauspec(ipatch,iwave) = tau1 +((ptau - p1) * &
                        dtauc(ilayer) / patch(ipatch)%atm(ilayer)%dp)
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
       
       spectrum = spectrum + (upflux*patch(ipatch)%cover)

    end do ! patch loop

       
    deallocate(upflux)


    
  end subroutine run_RT

    
  
end module setup_RT
