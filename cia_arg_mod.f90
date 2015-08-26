module cia_arg


contains



  subroutine get_cia(cia,ciatemp,grav,ch4index)

    use sizes
    use common_arrays
    use define_types
    use phys_const
    use atmos_ops
    
    implicit none
    real,dimension(nciatemps,nwave)::ph2h2,ph2He,ph2h,ph2ch4
    real,dimension(4,nciatemps,nwave),intent(in):: cia
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
       ! number density in amagats P0 in millibar
       n_amg = (patch(1)%atm(ilayer)%press / 1013.25 ) * (273.15 / patch(1)%atm(ilayer)%temp)

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

  
end module cia_arg
