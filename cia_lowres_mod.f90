module cia_lowres


contains



  subroutine get_cia_LR(grav,ch4index)

    use sizes
    use common_arrays
    use define_types
    use phys_const
    use atmos_ops
    
    implicit none
    double precision,dimension(nciatemps,ncwave)::ph2h2,ph2He,ph2h,ph2ch4
    double precision, dimension(ncwave) :: logciah2h2_LR, logciah2he_LR,logciah2ch4_LR   
    double precision, dimension(nwave) :: ciah2h2, ciah2he,ciah2ch4    
    real, dimension(nciatemps) :: ciatemp,tdiff
    double precision,dimension(ncwave)::ciawaven, wdiff
    real :: junk, grav
    double precision :: intfact, n_amg
    integer :: ilayer, idum1, idum2, ch4index,tcia1,tcia2, icwaven,iciatemp,iwave,oldw1,oldw2
    real :: fboth, fH2, fHe
    character(len=50):: ciafile1

    ! This code calculates the optical depth due to CIA from H2-H2, H2-He and H2-CH4, neglects H2-H
    ! It also neglects the other bound-free opacity we might expect
    ! This version interpolates in low-res version of CIA table, then rebins to working resolution at end
    ! before optical depth is calculated


      
    write(ciafile1,"(A,I0)") "../LineLists/final1_abel_CIA.dat"

          
    open(15,file=ciafile1,status="old")
    read(15,*) idum1, idum2

    if (idum1 .ne. ncwave .or. idum2 .ne. nciatemps) then
       write(*,*) " Problem with low-res CIA table : ", trim(ciafile1)
       stop
    end if
    
    do iciatemp = 1, nciatemps
       read(15,*) ciatemp(iciatemp)
       do icwaven = 1, ncwave
          read(15,*) ciawaven(icwaven), ph2h2(iciatemp,icwaven), ph2he(iciatemp,icwaven), &
               ph2h(iciatemp,icwaven), ph2ch4(iciatemp,icwaven)
       end do
    end do
    close(15)
    
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


          
          logciaH2H2_LR = (((ph2h2(2,:) - ph2h2(1,:))* intfact) + ph2h2(1,:))
          logciaH2He_LR = (((ph2he(2,:) - ph2he(1,:))* intfact) + ph2he(1,:))
          logciaH2CH4_LR = (((ph2ch4(2,:) - ph2ch4(1,:))* intfact) + ph2ch4(1,:))
       else
          intfact =  (log10(patch(1)%atm(ilayer)%temp) - log10(ciatemp(tcia1))) / (log10(ciatemp(tcia2)) - log10(ciatemp(tcia1)))
          
          ! resultant CIA is then:
          logciaH2H2_LR = (((ph2h2(tcia2,:) - ph2h2(tcia1,:))* intfact) + ph2h2(tcia1,:))
          logciaH2He_LR = (((ph2he(tcia2,:) - ph2he(tcia1,:))* intfact) + ph2he(tcia1,:))
          logciaH2CH4_LR = (((ph2ch4(tcia2,:) - ph2ch4(tcia1,:))* intfact) + ph2ch4(tcia1,:))


       endif

       ! Now we need to resample to higher resolution

       do iwave= 1 , nwave

          wdiff = abs(ciawaven - wavenum(iwave))

          oldw1 = minloc(wdiff,1)

          if (ciawaven(oldw1) .lt. wavenum(iwave)) then
             oldw2 = oldw1 + 1
          else
             oldw2 = oldw1
             oldw1 = oldw2 - 1
          end if
          
          
          intfact = (log10(wavenum(iwave)) - log10(ciawaven(oldw1))) / &
               (log10(ciawaven(oldw2)) - log10(ciawaven(oldw1)))

          ciaH2H2(iwave) = 10.0**(((logciaH2H2_LR(oldw2) - logciaH2H2_LR(oldw1))*intfact) + logciaH2H2_LR(oldw1))
          ciaH2He(iwave) = 10.0**(((logciaH2He_LR(oldw2) - logciaH2He_LR(oldw1))*intfact) + logciaH2He_LR(oldw1))
          ciaH2CH4(iwave) = 10.0**(((logciaH2CH4_LR(oldw2) - logciaH2CH4_LR(oldw1))*intfact) + logciaH2CH4_LR(oldw1))

       end do ! wave do

       
       ! Neglecting H2-H CIA  
       ! number density in amagats P0 in millibar
       n_amg = (patch(1)%atm(ilayer)%press / 1013.25 ) * (273.15/patch(1)%atm(ilayer)%temp)
       
       ! now calculate optical depth using amagats density unit
       patch(1)%atm(ilayer)%opd_cia = (n_amg**2 * patch(1)%atm(ilayer)%fH2 * patch(1)%atm(ilayer)%dz) * &
            ((patch(1)%atm(ilayer)%fH2 * ciaH2H2) &
            + (patch(1)%atm(ilayer)%fHe * ciaH2He) &
            + (patch(1)%atm(ilayer)%gas(ch4index)%VMR * ciaH2CH4))
       
    end do ! layer do
    
  end subroutine get_cia_LR

  
end module cia_lowres
