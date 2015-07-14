module setup_twostr

contains

  subroutine run_twostr(spectrum,nw1,nw2)

    use sizes
    use common_arrays
    use define_types
    use phys_const
    use atmos_ops

    
    implicit none

    integer, parameter :: nlevel = nlayers+1
    integer, parameter :: MAXCLY = nlayers
    integer, parameter :: MAXULV = 1


    ! now set up where we want the fluxes.
    ! we just want the top of atmosphere to start
    ! can worry about bonus bits for tracking contributions later
    
    integer ::ntau
    double precision,dimension(MAXULV) :: utau
    
    double precision,dimension(nwave), INTENT(OUT):: spectrum
    double precision, dimension(nlayers) :: DTAUC, SSALB, GG
    double precision, dimension(nlayers+1) :: temper
    double precision :: WVNMLO, WVNMHI, wint
    integer :: ipatch,ilayer, iwave
    double precision:: umu0
    double precision,dimension(MAXULV) :: RFLDIR,RFLDN,FLUP,DFDT,UAVG
    double precision:: BTEMP, TTEMP
    double precision,dimension(nwave):: upflux
    logical, dimension(2) :: PRNT
    character(len=127):: HEADER
    integer :: NLYR, IBCND, ierror
    logical :: PLANK,  USRTAU, spher, quiet, deltam
    double precision :: FBEAM, FISOT,  ALBEDO , TEMIS
    double precision:: radius, zd
 

    integer :: nw1, nw2
    
    HEADER = ""

    PRNT = [.FALSE., .FALSE.]

    QUIET = .false.

    deltam = .true.
    ! set values for non-parameter disort input
    NLYR = nlayers
    NTAU = 1
    UTAU = 1e-5
    radius = 0.
    spher = .FALSE.
    zd = 0.
    
    USRTAU    = .TRUE.  ! return quantities at each pre-defined layer
    IBCND     = 0        ! boundary conditions, change if Planck is false
    FBEAM     = 0.0      ! parallel beam at the top. 0 for stars
    ! if this is not 0, need to specify umu0 and phi0
    FISOT     = 0.0      ! not sure, was 1.0 / PI
    ALBEDO = 0.0001E0
    PLANK    = .TRUE.    ! need this to use temperature structure
    TEMIS = 0.01         ! need to give top layer a bit of emissivity 
    

    ! all temps the same across patches
    call set_temp_levels(patch(1)%atm%temp,temper)

    spectrum = 0.0
    
    do ipatch = 1, npatch
       ! add up the taus to get extinction
       do ilayer = 1, nlayers
        
          patch(ipatch)%atm(ilayer)%opd_ext = &
               patch(ipatch)%atm(ilayer)%opd_ext + &
               patch(ipatch)%atm(ilayer)%opd_lines + &
               patch(ipatch)%atm(ilayer)%opd_CIA
       end do
     



       upflux = 0.0
       do iwave = nw2, nw1
       
          DTAUC = patch(ipatch)%atm%opd_ext(iwave)              
          
          if (patch(ipatch)%cloudy) then
             GG = patch(ipatch)%atm%gg(iwave)       
             SSALB = patch(ipatch)%atm%opd_scat(iwave) / &
                  patch(ipatch)%atm%opd_ext(iwave)

          else
             GG = 1.
             SSALB = 0.0
          end if
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
          
          write(*,*) "HERE!"
          
          call twostr( albedo, btemp, deltam, dtauc, fbeam, fisot,& 
               gg, header, ierror, maxcly, maxulv, nlyr, plank,&
               ntau, prnt, quiet, radius, spher, ssalb, temis,&  
               temper, ttemp, umu0,  usrtau, utau, wvnmlo,&
               wvnmhi, zd, dfdt, flup, rfldir, rfldn, uavg )


         
          
          ! convert to flux density W/m2/um
          ! need interval in um not cm^-1
          wint = (1.0e4 / wavenum(iwave)) * ((WVNMHI - WVNMLO)/ wavenum(iwave))
          
          upflux(iwave) = FLUP(1) / wint

          write(*,*) "twostr done. IERROR = ",ierror
          
       end do ! wave loop
       
       spectrum = spectrum + (upflux*patch(ipatch)%cover)

       
    end do ! patch loop

  end subroutine run_twostr

    
  
end module setup_twostr
