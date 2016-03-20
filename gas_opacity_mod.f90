module gas_opacity
  

contains

  subroutine line_mixer(linelist)

    use sizes
    use common_arrays
    use phys_const
    use define_types


    implicit none

    
    double precision,intent(in):: linelist(:,:,:,:)
    integer :: Tlay1, Tlay2, torder, ilayer, iounit
    double precision, allocatable,dimension(:,:) :: logkap1, logkap2,logintkappa
    real, dimension(nlinetemps):: tdiff
    double precision :: ndens, intfact, junk
    character(len=50) :: lines1,lines2,name
    
    ! counters
    integer:: igas, j
    
    allocate(logkap1(ngas,nwave),logkap2(ngas,nwave),logintkappa(ngas,nwave))
    
    ! get line temp array locations bracketing our temperature
    
    do ilayer= 1, nlayers
       
       tdiff = abs(linetemps - patch(1)%atm(ilayer)%temp)
    
    
       Tlay1 = minloc(tdiff,1)
    
    
       if (linetemps(Tlay1) .lt. patch(1)%atm(ilayer)%temp) then
          Tlay2 = Tlay1+1
       else
          Tlay2 = Tlay1 - 1
       end if
       
       
       
       
       if (patch(1)%atm(ilayer)%temp .lt. linetemps(1)) then
          Tlay1 = 1
          Tlay2 = 2
       else if (patch(1)%atm(ilayer)%temp .gt. linetemps(nlinetemps)) then
          Tlay1 = nlinetemps - 1
          Tlay2 = nlinetemps
       endif
       
       
       
       ! get linear interpolation factor
       
       if (Tlay1 .gt. Tlay2) then
          torder = 1
          intfact = (log10(patch(1)%atm(ilayer)%temp) - &
               log10(linetemps(Tlay2))) / &
               (log10(linetemps(Tlay1)) - log10(linetemps(Tlay2)))
       else
          torder = 2
          intfact =  (log10(patch(1)%atm(ilayer)%temp) -&
               log10(linetemps(Tlay1))) / &
               (log10(linetemps(Tlay2)) - log10(linetemps(Tlay1)))
       endif
       
       if (patch(1)%atm(ilayer)%temp .gt. linetemps(nlinetemps)) then
          intfact =  (log10(patch(1)%atm(ilayer)%temp) &
               - log10(linetemps(Tlay2))) &
               / (log10(linetemps(Tlay2)) - log10(linetemps(Tlay1)))
       endif
       
       
       
       patch(1)%atm(ilayer)%opd_lines = 0.0
       
       
       !! CODE FUCKS UP HERE.. segfault 
       logkap1 = linelist(:,ilayer,Tlay1,:)
       logkap2 = linelist(:,ilayer,Tlay2,:) 
       
       !  do the interpolation in log(cross section) !!!         !


       if (torder .eq. 1) then 
          logintkappa = ((logkap1 - logkap2)*intfact)+logkap2
       else if (torder .eq. 2) then
          logintkappa = ((logkap2 - logkap1)*intfact)+logkap1
       else
          write(*,*) "something wrong with interpolate order"           
       endif
       
       if (patch(1)%atm(ilayer)%temp .gt. linetemps(nlinetemps)) then
          logintkappa = ((logkap2 - logkap1)*intfact)+logkap2
       endif

       ! now we've got cross section for gas in layer - intkappa
       ! (cm2 / molecule)
       ! we want the optical depth.. 
       ! so we multiply by number density
       !(which is in /m3 to get column density ( / cm)
       ! and sum over the layer thickness (which we calculated in METRES)
       ! to get optical depth

       
       do igas = 1, ngas       
          patch(1)%atm(ilayer)%opd_lines = patch(1)%atm(ilayer)%opd_lines + &
               ((patch(1)%atm(ilayer)%gas(igas)%VMR * &
               patch(1)%atm(ilayer)%ndens * &
               patch(1)%atm(ilayer)%dz * 1d-4) &
               * (10.0**logintkappa(igas,:)))
       end do
    end do ! do loop
 
    deallocate(logkap1, logkap2,logintkappa)
    
  end subroutine line_mixer




  subroutine get_ray(ch4index)
    
    use sizes
    use common_arrays
    use phys_const
    use define_types
    
    
    implicit none
    
    integer :: ilayer,ch4index,nn,iwave,ng
    double precision :: cfray, XN0,cold
    double precision, dimension(3):: gasss,dpol
    double precision,dimension(2,3):: gnu
    double precision,allocatable,dimension(:):: wa,tec, taur
    dpol = [1.022d0,1.0d0,1.0d0]
    
    allocate(wa(nwave),tec(nwave),taur(nwave))
    
    gnu(1,1)=1.355d-4
    gnu(2,1)=1.235d-6
    gnu(1,2)=3.469d-5
    gnu(2,2)=8.139d-8
    gnu(1,3)=4.318d-4
    gnu(2,3)=3.408d-6
    XN0=2.687d19
    
    cfray = 32.d0*pi**3*1.d21/(3.d0*2.687d19) 
    wa = wavelen !* 1e-4
   
    do ilayer = 1, nlayers    
       gasss(1) = patch(1)%atm(ilayer)%fH2
       gasss(2) = patch(1)%atm(ilayer)%fHe
       if (ch4index .ne. 0) then
          gasss(3) = patch(1)%atm(ilayer)%gas(ch4index)%VMR
          ng = 3
       else
          ng = 2
       end if
       
       
       
       cold = patch(1)%atm(ilayer)%ndens * patch(1)%atm(ilayer)%dz * 1.e-4
    
    
       taur = 0.0      
       do nn =1,ng
          tec = cfray*(dpol(nn)/wa**4)*(gnu(1,nn)+gnu(2,nn) &
               / wa**2)**2
          taur = taur + COLD*gasss(nn) * tec * 1.d-5 / XN0
          
       end do !nn loop

       !    taur = cold * 8.49e-45 / wa**4
       patch(1)%atm(ilayer)%opd_rayl = taur
    end do ! layer do
    
    deallocate(tec,taur,wa)
    
  end subroutine get_ray



  subroutine get_cia(cia,ciatemp,grav,ch4index)

    use sizes
    use common_arrays
    use define_types
    use phys_const
    use atmos_ops
    
    implicit none
    real,dimension(nciatemps,nwave)::ph2h2,ph2He,ph2h,ph2ch4
    real,intent(inout):: cia(:,:,:)
    real, dimension(nwave) :: ciah2h2, ciah2he,ciah2ch4    
    real,dimension(nciatemps), intent(in) :: ciatemp 
    real, intent(in) :: grav
    integer, intent(in) :: ch4index
    real,dimension(nciatemps) :: tdiff
    real :: intfact, n_amg
    integer :: ilayer, idum1, idum2,tcia1,tcia2, icwaven,iciatemp
    real :: fboth, fH2, fHe

    ! This code calculates the optical depth due to CIA from H2-H2, H2-He and H2-CH4, neglects H2-H
    ! It also neglects the other bound-free opacity we might expect
    ! This version interpolates in low-res version of CIA table, then rebins to working resolution at end
    ! before optical depth is calculated


    ph2h2  = cia(1,:,:)
    ph2he  = cia(2,:,:)
    ph2h   = cia(3,:,:)
    ph2ch4 = cia(4,:,:)


    do ilayer = 1, nlayers
       

       ! now we need to inpterpolate for the temperature


       
       tdiff = abs(ciatemp - patch(1)%atm(ilayer)%temp)
       
       tcia1 = minloc(tdiff,1)
       if (ciatemp(tcia1) .lt. patch(1)%atm(ilayer)%temp)  then
          tcia2 = tcia1 + 1
       else
          tcia2 = tcia1
          tcia1 = tcia2 - 1
       end if
       
       
       if (tcia1 .eq. 0) then
          intfact = (log10(patch(1)%atm(ilayer)%temp) - log10(ciatemp(1))) / (log10(ciatemp(2)) - log10(ciatemp(1)))
          
          ciaH2H2 = 10.0**(((ph2h2(2,:) - ph2h2(1,:))* intfact) + ph2h2(1,:))
          ciaH2He = 10.0**(((ph2he(2,:) - ph2he(1,:))* intfact) + ph2he(1,:))
          ciaH2CH4 = 10.0**(((ph2ch4(2,:) - ph2ch4(1,:))* intfact) + ph2ch4(1,:))
       else
          intfact =  (log10(patch(1)%atm(ilayer)%temp) - log10(ciatemp(tcia1))) / (log10(ciatemp(tcia2)) - log10(ciatemp(tcia1)))
          
          ! resultant CIA is then:
          ciaH2H2 = 10.0**(((ph2h2(tcia2,:) - ph2h2(tcia1,:))* intfact) + ph2h2(tcia1,:))
          ciaH2He = 10.0**(((ph2he(tcia2,:) - ph2he(tcia1,:))* intfact) + ph2he(tcia1,:))
          ciaH2CH4 = 10.0**(((ph2ch4(tcia2,:) - ph2ch4(tcia1,:))* intfact) + ph2ch4(tcia1,:))
          
          
       endif
       
       
       ! Neglecting H2-H CIA  
       ! number density in amagats P0 in bar
       n_amg = (patch(1)%atm(ilayer)%press / 1.01325) * (273.15 / patch(1)%atm(ilayer)%temp)
       
       ! now calculate optical depth using amagats density unit
       ! put layer thickness in cm as CIA cross sections are in cm^-1 amg^-2
       
       ! Check if CH4 is present and include CH4-H2 if present       
       if (ch4index .ne. 0) then 
          patch(1)%atm(ilayer)%opd_cia = (n_amg**2. * patch(1)%atm(ilayer)%fH2 * patch(1)%atm(ilayer)%dz*100.) * &
               ((patch(1)%atm(ilayer)%fH2 * ciaH2H2) &
               + (patch(1)%atm(ilayer)%fHe * ciaH2He) &
               + (patch(1)%atm(ilayer)%gas(ch4index)%VMR * ciaH2CH4))
       else
          ! If no CH4 just sum the H2-H2 and H2-He CIAs
          patch(1)%atm(ilayer)%opd_cia = (n_amg**2. * patch(1)%atm(ilayer)%fH2 * patch(1)%atm(ilayer)%dz*100.) * &
               ((patch(1)%atm(ilayer)%fH2 * ciaH2H2) &
               + (patch(1)%atm(ilayer)%fHe * ciaH2He))
       end if
       
    end do ! layer do
  end subroutine get_cia
  
  
end module gas_opacity
