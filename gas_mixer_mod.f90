module gas_mixing

  implicit none


contains

  subroutine line_mixer(layer,mu,grav,opd_lines)

    use sizes
    use phys_const
    use define_types

    real,intent(IN) :: grav, mu
    type(a_layer), intent(IN) :: layer
    double precision, intent(OUT) :: opd_lines
    integer :: Tlay1, Tlay2, torder
    double precision, dimension(nwave) :: kappa1,kappa2,intkappa,totkappa
    double precision, dimension(nwave) :: logkap1, logkap2
    real, dimension(nlinetemps):: linetemp,tdiff
    double precision :: ndens
    character(len=50) :: lines1,lines2
    
    
    ! set up the line temp array
    
    call set_line_temps

    ! get line temp array locations bracketing our temperature

    tdiff = abs(linetemp - layer%temp)

    
    Tlay1 = minloc(tdiff,1)

      
    if (linetemp(Tlay1) .lt. T) then
       Tlay2 = Tlay1+1
    else
       Tlay2 = Tlay1 - 1
    endif


      ! get linear interpolation factor

      if (Tlay1 .gt. Tlay2) then
         torder = 1
         intfact = (layer%temp - linetemp(Tlay2)) / (linetemp(Tlay1) - linetemp(Tlay2))
         ! TK test line
         write(*,*) "torder = ",torder
         write(*,*) linetemp(Tlay2),linetemp(Tlay1)
      else
         torder = 2
         infact =  (layer%temp -linetemp(Tlay1)) / (linetemp(Tlay2) - linetemp(Tlay1))
         ! TK test line
         write(*,*) "torder = ",torder
         write(*,*) linetemp(Tlay2),linetemp(Tlay1)
      endif


      ! now get the files by their Tlay and pressure locations
      ! gases are identified by name (character!)
      ! this will allow the relevant linelists to be grabbed by T and P index



      ! loope through the gases

      totkappa = 0.0

      do i = 1, ngas

         write(lines1,"(A,A,A1,A,A1,I0,A1,I0)")"../LineLists/",trim(layer%gas(i)%name),"/",trim(layer%gas(i)%name),"_",Tlay1,"_",layer%index
         write(lines2,"(A,A,A1,A,A1,I0,A1,I0)") "../LineLists/",trim(layer%gas(i)%name),"/",trim(layer%gas(i)%name),"_",Tlay2,"_",layer%index

         ! TK test line
         write(*,*) lines1
         write(*,*) lines2



       ! now read in the line lists and sum, weighted by abudance/fraction

         
         open(15,file=lines1)
         do j =1, listheadlines
            read(15,*)
         enddo
         do j = 1,nwave
            read(15,*) wavenum(j),kappa1(j)
         enddo
         close(15)
         open(15,file=lines2)
         do j =1, listheadlines
            read(15,*)
         enddo
         do j = 1,nwave
            read(15,*) junk,kappa2(j)
         enddo
         close(15)

         !  do the interpolation in log(cross section) !!!         !
         logkap1 = log10(kappa1)
         logkap2 = log10(kappa2)
         if (torder .eq. 1) then 
            intkappa = ((logkap1 - logkap2)*intfact)+logkap2
         else if (torder .eq. 2) then
            intkappa = ((logkap2 - logkap1)*intfact)+logkap1
         else
            write(*,*) "something wrong with interpolate order"           
         endif
         totkappa = totkappa + (layer%gas(i)%VMR * 10**intkappa)
      enddo

      ! now we've got total cross section for layer - totkappa (cm2 / molecule)
      ! we want the optical depth.. 
      ! so we multiply by number density to get column density ( / cm)
      ! and sum over the layer thickness (which we calculated in METRES

      ! number density in /m3:

      ndens = layer%press  / (K_BOLTZ * layer%temp)
      
      ! optical depth
      
      opd_lines = totkappa * 10**(-4) * ndens * layer%dz 


      

  end subroutine line_mixer


end module gas_mixing
