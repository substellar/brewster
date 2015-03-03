program opacity_mixer

      implicit double precision (a-h,o-z)

      !include 'globals.h'

      integer:: ngas, Tlay1, Tlay2, Play,torder
      integer, allocatable :: gas(:)
      real, allocatable::abgas(:)
      real :: T ,junk,intfact
      real, dimension(319188) :: wavenum,lambda
      double precision, dimension(319188) :: kappa1,kappa2,intkappa,totkappa
      character(len=50) :: abfile,lines1,lines2
      ! temp and pressure arrays
      real, dimension(27) :: temp,tdiff
      real, dimension(16) :: pres

!     This code reads in line lists and mixes them for a layer
!     Pressure layers are set to same as in supplied line lists for speed
!     T is interpolated

      nwave = 319188
      write(*,*) "give number of gases"
      read(*,*) ngas
      write(*,*) "Give gas list file"
      read(*,*) abfile

      write(*,*) "Give T and Player (latter is array location)"
      read(*,*) T, Play

      
      allocate (gas(ngas),abgas(ngas))

      open(15,file=abfile,status='old')

      do i =1,ngas
         read(15,*) gas(i), abgas(i)
      enddo

      close(15)

      

      ! Now we want to find the array locations for our temperature point
      ! before setting up a loop that will read in linelists for each gas
      ! and get list for layer

      ! temp array
      temp=[75., 85., 100., 120., 140., 160., 180., 200., 230., 260., &
           300., 350., 400.,500., 650., 800., 1000., 1200., 1400., 1600.,&
           1800., 2000., 2300., 2600.,3000., 3500., 4000.]

      ! pressures

      pres=[1.e-02, 3.e-02, 1.e-1, 3.e-1,1.e0, 3.e0,1.e1, 3.e1, 1.e2, &
           3.e2,1.e3, 3.e3,1.e4, 3.e4,1.e5,3.e5] 


      ! get Temp locations

      tdiff = abs(temp - T)
      write(*,*) tdiff

      Tlay1 = minloc(tdiff,1)

      
      if (temp(Tlay1) .lt. T) then
         Tlay2 = Tlay1+1
      else
         Tlay2 = Tlay1 - 1
      endif


      ! get linear interpolation factor

      if (Tlay1 .gt. Tlay2) then
         torder = 1
         intfact = (T - temp(Tlay2)) / (temp(Tlay1) - temp(Tlay2))
         write(*,*) "torder = ",torder
         write(*,*) temp(Tlay2),temp(Tlay1)
      else
         torder = 2
         infact =  (T -temp(Tlay1)) / (temp(Tlay2) - temp(Tlay1))
         write(*,*) "torder = ",torder
         write(*,*) temp(Tlay2),temp(Tlay1)
      endif

      
      ! now get the files by their Tlay, Play locations
      ! gases will be identified by integer gas
      ! this will allow the relevant linelists to be grabbed by T and P index


      ! gases are:

      ! 1: H2O
      ! 2: CH4
      ! more to come
      totkappa = 0.0
      
      do i = 1, ngas
         
         if (gas(i) .eq. 1) then
            write(lines1,"(A11,I0,A1,I0)") "../H2O/h2o_",Tlay1,"_",Play
            write(lines2,"(A11,I0,A1,I0)") "../H2O/h2o_",Tlay2,"_",Play
         else if (gas(i) .eq. 2) then
            write(lines1,"(A11,I0,A1,I0)") "../CH4/ch4_",Tlay1,"_",Play
            write(lines2,"(A11,I0,A1,I0)") "../CH4/ch4_",Tlay2,"_",Play
         endif

!         write(*,*) lines1
!         write(*,*) lines2
         
       ! now read in the line lists and sum, weighted by abudance/fraction


         open(15,file=lines1)
         do j =1, 23
            read(15,*)
         enddo
         do j = 1,nwave
            read(15,*) wavenum(j),kappa1(j)
         enddo
         close(15)
         open(15,file=lines2)
         do j =1, 23
            read(15,*)
         enddo
         do j = 1,nwave
            read(15,*) junk,kappa2(j)
         enddo
         close(15)
         if (torder .eq. 1) then 
            intkappa = ((kappa1 - kappa2)*intfact)+kappa2
         else if (torder .eq. 2) then
            intkappa = ((kappa2 - kappa1)*intfact)+kappa1
         else
            write(*,*) "something wrong with interpolate order"           
         endif
         totkappa = totkappa + (abgas(i)*intkappa)
      enddo

      open(16,file="testopacity.dat",status='new')
      write(16,*) "wavenum mixedkappa"
      do i = 1, nwave
         write(16,*) wavenum(i),totkappa(i)
      enddo
      close(16)
      
      end program
