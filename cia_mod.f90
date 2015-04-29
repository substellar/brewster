module cia


contains



  subroutine get_cia(grav,ch4index)

    use sizes
    use common_arrays
    use define_types
    use phys_const
    use atmos_ops
    
    implicit none

    real:: p1, p2, temp1,temp2
    real, dimension(nlayers+1):: Tlevels(0:nlayers) 
    double precision,dimension(nwave)::ph2h2_1,ph2He_1,ph2h_1,ph2ch4_1,ph2h2_2,ph2He_2,ph2h_2,ph2ch4_2
    double precision, dimension(nwave) :: ciah2h2, ciah2he,ciah2ch4
    real, dimension(nciatemps) :: ciatemp,tdiff
    real :: junk, grav, intfact, acoef,bcoef,coef1,part1,part2,part3,part4
    integer :: i, j, idum2, ch4index,tcia1,tcia2
    real :: fboth, fH2, fHe
    character(len=50):: ciafile1,ciafile2
    ! we'll need level temps

    call set_temp_levels(atm%temp,Tlevels)

    open(10,file="../LineLists/CIA/cia_temps.dat",status="old")
    read(10,*) ciatemp
      
     
    do i = 1, nlayers
       
       temp1 = Tlevels(i-1)
       temp2 = Tlevels(i)
       ! get the pressure levels at boundaries:
       
       if (i .eq. nlayers) then
          p1 =  exp((0.5)*(log(atm(i-1)%press * atm(i)%press)))
          p2 = atm(i)%press**2 / p1
       else
          p1 = exp(((1.5)*log(atm(i)%press)) - ((0.5)*log(atm(i+1)%press)))
          p2 = exp((0.5)*(log(atm(i)%press * atm(i+1)%press)))
       endif
       
       ! first thing is coefficient calculation lifted from Mark's code:
       
       ! put p's into bar (they're in millibar to start)
       
       p1 = p1/1000.0
       p2 = p2/1000.0
       
       ! SET UP THE COEFFICIENT TO REDUCE MASS PATH TO STP ...SEE NOTES
       ! T0 =273.15   PO=1.01325 BAR
       ! NOTE THAT XMU IS ALREADY A LAYER AVERAGE QUANTITY
       !    McKay 1/96 derivation:
       part1= (atm(i)%temp/(temp1*temp2))
       part2=(temp2*p2 - temp1*p1)
       part3=(p2-p1)
       part4=part1*part2/part3
       ACOEF = (atm(i)%temp/(temp1*temp2))*(temp2*p2 - temp1*p1)/(p2-p1)
       
       BCOEF = (atm(i)%temp/(temp1*temp2))* (temp1 - temp2)/(p2-p1)

       ! convert RGAS to units expected here...ergs/mol/K.
       ! also check gravity and mu units
       ! also gravity needs to be in cm/s
       
       COEF1=(R_GAS)*273.15**2*.5E5* (ACOEF* (p2**2 - p1**2) + &
            BCOEF*(2./3.)*(p2**3 - p1**3) ) &
            /(1.01325**2 *(100.0*grav)*atm(i)%temp*atm(i)%mu)


       ! now we need to inpterpolate for the temperature

       tdiff = abs(ciatemp - atm(i)%temp)

       tcia1 = minloc(tdiff,1)

       if (ciatemp(tcia1) .lt. atm(i)%temp)  then
          tcia2 = tcia1 + 1
       else
          tcia2 = tcia1
          tcia1 = tcia2 - 1
       end if

       if (tcia1 .eq. 0) then
          intfact = (log10(atm(i)%temp) - log10(ciatemp(1))) / (log10(ciatemp(2)) - log10(ciatemp(1)))

          write(ciafile1,"(A,I0)") "../LineLists/CIA/cia_highres_tempindex_1"
          write(ciafile2,"(A,I0)") "../LineLists/CIA/cia_highres_tempindex_2"
          
          open(15,file=ciafile1,status="old")
          read(15,*) junk, idum2
          
          if (junk .ne. ciatemp(1)) then
             write(*,*) "CIA file input error: ", trim(ciafile1)
             stop
          end if
          
          do j = 1, nwave
             read(15,*) junk, ph2h2_1(j),ph2he_1(j),ph2h_1(j),ph2ch4_1(j)
          end do
          close(15)
          
          open(15,file=ciafile2,status="old")
          read(15,*) junk, idum2
          if (junk .ne. ciatemp(2)) then
             write(*,*) "CIA file input error: ", trim(ciafile2)
             stop
          end if
          do j = 1, nwave
             read(15,*) junk, ph2h2_2(j),ph2he_2(j),ph2h_2(j),ph2ch4_2(j)
          end do
          close(15)
          ciaH2H2 = 10.0**(((ph2h2_2 - ph2h2_1)* intfact) + ph2h2_1)
          ciaH2He = 10.0**(((ph2he_2 - ph2he_1)* intfact) + ph2he_1)
          ciaH2CH4 = 10.0**(((ph2ch4_2 - ph2ch4_1)* intfact) + ph2ch4_1)
       else
          intfact =  (log10(atm(i)%temp) - log10(ciatemp(tcia1))) / (log10(ciatemp(tcia2)) - log10(ciatemp(tcia1)))
          
          ! resultant CIA is then:
          write(ciafile1,"(A,I0)") "../LineLists/CIA/cia_highres_tempindex_",tcia1
          write(ciafile2,"(A,I0)") "../LineLists/CIA/cia_highres_tempindex_",tcia2
          
          open(15,file=ciafile1,status="old")
          read(15,*) junk, idum2
          
          if (abs(junk - ciatemp(tcia1)) .gt. 10.0) then
             write(*,*) "CIA file input error: ", trim(ciafile1)
             stop
          end if
          
          do j = 1, nwave
             read(15,*) junk, ph2h2_1(j),ph2he_1(j),ph2h_1(j),ph2ch4_1(j)
          end do
          close(15)
          
          open(15,file=ciafile2,status="old")
          read(15,*) junk, idum2
          if (junk .ne. ciatemp(tcia2)) then
             write(*,*) "CIA file input error: ", trim(ciafile2)
             stop
          end if
          do j = 1, nwave
             read(15,*) junk, ph2h2_2(j),ph2he_2(j),ph2h_2(j),ph2ch4_2(j)
          end do
          close(15)
          ciaH2H2 = 10.0**(((ph2h2_2 - ph2h2_1)* intfact) + ph2h2_1)
          ciaH2He = 10.0**(((ph2he_2 - ph2he_1)* intfact) + ph2he_1)
          ciaH2CH4 = 10.0**(((ph2ch4_2 - ph2ch4_1)* intfact) + ph2ch4_1)


       endif
       ! Neglecting H2-H CIA  
       
       atm(i)%opd_cia = COEF1 * ((atm(i)%fH2 * atm(i)%fH2 * ciaH2H2) + &
            (atm(i)%fH2 * atm(i)%fHe * ciaH2He) + &
            (atm(i)%fH2 * atm(i)%gas(ch4index)%VMR * ciaH2CH4))
       
       
       
    end do
    
  end subroutine get_cia

  
end module cia
