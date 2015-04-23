module setup_disort




contains

  subroutine run_disort(spectrum,do_clouds,nw1,nw2)

    use sizes
    use common_arrays
    use define_types
    use phys_const
    use atmos_ops
!    use disort_params    
!    use disort_mod
    
    implicit none

    integer, parameter :: nlevel = nlayers+1
    integer, parameter :: MAXCLY = nlayers
    integer, parameter :: MAXMOM = 8
    integer, parameter :: MAXPHI = 3
    integer, parameter :: MAXULV = 1
    integer, parameter :: MAXUMU = 10


    ! now set up where we want the fluxes.
    ! we just want the top of atmosphere to start
    ! can worry about bonus bits for tracking contributions later
    
    integer ::ntau
    real,dimension(MAXULV) :: utau
    
    
    logical, intent(IN) :: do_clouds
    real,dimension(nwave), INTENT(OUT):: spectrum
    real, dimension(nlayers) :: DTAUC, SSALB
    real, dimension(nlayers+1) :: temper
    real :: WVNMLO, WVNMHI
    integer :: i,j
    real,dimension(0:MAXMOM,nlayers) :: PMOM
    real:: phi(maxphi),phi0, umu(maxumu), umu0
    real,dimension(MAXULV) :: RFLDIR,RFLDN,FLUP,DFDT,UAVG
    real,dimension(maxumu,maxulv,maxphi) :: UU
    real,dimension(maxumu) :: ALBMED, TRNMED
    real:: BTEMP, TTEMP
    real,dimension(nwave):: upflux
    logical, dimension(5) :: PRNT
    character(len=127):: HEADER
    integer :: NUMU,NSTR,NMOM,NLYR, NPHI, IBCND
    logical :: LAMBER,  PLANK,  USRTAU, USRANG, ONLYFL
    real :: FBEAM, FISOT,  ALBEDO , ACCUR, TEMIS



    


    
    integer :: nw1, nw2
    
    HEADER = ''

     PRNT = [.TRUE., .FALSE.,.FALSE.,.FALSE., .TRUE. ]

     ! set values for non-parameter disort input
     NSTR = 8
     NMOM = 8
     NUMU = 8 ! same as NSTR
     NLYR = nlayers
     NPHI      = 0
     NTAU = 1
     UTAU = 1e-5

     USRTAU    = .TRUE.  ! return quantites at each pre-defined layer
     USRANG    = .FALSE.
     ONLYFL    = .TRUE.   ! return only fluxes      
     IBCND     = 0        ! boundary conditions, change if Planck is false
     FBEAM     = 0.0      ! parallel beam at the top. 0 for stars
     ! if this is not 0, need to specify umu0 and phi0
     FISOT     = 0.0      ! not sure, was 1.0 / PI
     LAMBER    = .TRUE.   ! don't care about lower boundary
     ALBEDO = 0.01E0
     ACCUR = 0.0 
     PLANK    = .TRUE.    ! need this to use temperature structure
     TEMIS = 0.01         ! need to give top layer a bit of emissivity 

     
     call set_temp_levels(atm%temp,temper)

     upflux = 0.0
  
    do i = nw2, nw1
       
       SSALB = atm%opd_scat(i) / atm%opd_ext(i)
       ! junk line to set SSALB in absence of dust
       SSALB = 0.0
       DTAUC = atm%opd_ext(i)

       ! set up wavenumber interval....

      
       if (i .eq. 1) then 
          WVNMLO = wavenum(1) - 0.5*(wavenum(2) - wavenum(1))
          WVNMHI =  wavenum(1) + 0.5*(wavenum(2) - wavenum(1))
       else if (i .eq. nwave) then
          WVNMLO = wavenum(nwave) - 0.5*(wavenum(nwave) - wavenum(nwave-1))
          WVNMHI =  wavenum(nwave) + 0.5*(wavenum(nwave) - wavenum(nwave-1))
       else
          WVNMLO = wavenum(i) - 0.5*(wavenum(i) - wavenum(i-1))
          WVNMHI = wavenum(i) + 0.5*(wavenum(i+1) -wavenum(i))
       end if

       ! set top and bottom boundary temperatures
       BTEMP = real(temper(nlayers+1))
       TTEMP = real(temper(1))


!       if (i .eq. nw2 + 3) STOP
       
       ! need PMOM
       ! get this calling Mark's GETMOM code
       ! 6 is Rayleigh+HG
       ! 2 ia just Rayleigh
       if (.not.do_clouds) then 
          DO j = 1, nlayers
             CALL GETMOM( 2, 0.0, NMOM, 0.0, PMOM(0:nmom,j) )
             ! fudge to fix the issue with PMOM values
          enddo
       else     
          do j = 1, nlayers
             ! TK test line
             write(*,*) "DOING GETMOM H-G"
             call GETMOM(6,atm(j)%gg(i),NMOM,PMOM(0,j))
          end do
       end if

 
       write(*,*) "Calling DISORT"
       call DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO, &
            WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, &
            UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, &
            FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS,&
            PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY,&
            MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, RFLDN,&
            FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )



       
       upflux(i) = FLUP(1)
       
    end do


    spectrum = upflux


    
  end subroutine run_disort

    
  
end module setup_disort
