module setup_disort



contains

  subroutine run_disort(spectrum,do_clouds,nw1,nw2)

    use sizes
    use common_arrays
    use define_types
    use phys_const
    use atmos_ops
    use disort_params    
    use disort_mod
    
    implicit none

    logical, intent(IN) :: do_clouds
    real,dimension(nwave), INTENT(OUT):: spectrum
    real, dimension(nlayers) :: DTAUC, SSALB
    real, dimension(nlayers+1) :: temper
    real :: WVNMLO, WVNMHI
    integer :: i,j
    real,dimension(0:NMOM,nlayers) :: PMOM
    real:: phi(nphi),phi0, umu(nstr), umu0
    real,dimension(ntau) :: RFLDIR,RFLDN,FLUP,DFDT,UAVG
    real,dimension(numu,ntau,nphi) :: UU
    real,dimension(numu) :: ALBMED, TRNMED
    real:: BTEMP, TTEMP
    real,dimension(nwave):: upflux
    logical, dimension(5) :: PRNT
    character(len=127):: HEADER

    integer :: nw1, nw2
    
    HEADER = ''

     PRNT = [.TRUE., .FALSE.,.FALSE.,.FALSE., .TRUE. ]

     call set_temp_levels(atm%temp,temper)

    upflux = 0.0
    ! TK test edit Just a subset of wavelength space... say 1 - 5um

    


     
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

       write(*,*) "Setup disort L75 BTEMP =", BTEMP
       write(*,*) "Setup disort L75 TTEMP =", TTEMP
       write(*,*) "Setup disort L75 temper = ", temper

       if (i .eq. nw2 + 1) STOP
       
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
!       call DISORT(nlayers, DTAUC, SSALB, NMOM, myPMOM, TEMPER, WVNMLO, &
!            WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, &
!            UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, &
!            FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, &
!            PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY,&
!       MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, RFLDN,&
!       FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )

       ! TK TEST LINE


       
       call DISORT(DTAUC, SSALB, PMOM, TEMPER, WVNMLO,&
            WVNMHI, UMU, PHI, UMU0, PHI0, BTEMP, TTEMP, &
            PRNT, HEADER, RFLDIR, RFLDN, &
            FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )

       
       upflux(i) = FLUP(1)
       
    end do


    spectrum = upflux


    
  end subroutine run_disort

    
  
end module setup_disort
