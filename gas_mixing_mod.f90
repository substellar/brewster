module gas_mixing



  implicit none


contains

  subroutine line_mixer(layer,opd_lines,index)

    use sizes
    use phys_const
    use define_types
    use common_arrays


    implicit none
    
    type(a_layer) :: layer
    double precision, dimension(nwave) :: opd_lines
    integer :: Tlay1, Tlay2, torder, index, iounit
    double precision, dimension(nwave) :: kappa1,kappa2,logintkappa,totkappa
    double precision, dimension(nwave) :: logkap1, logkap2
    real, dimension(nlinetemps):: tdiff
    double precision :: ndens, intfact, junk
    character(len=50) :: lines1,lines2,name

    ! counters
    integer:: i, j
    
    ! set up the line temp array - check this is working properly
    
    call set_line_temps
    
    
    ! get line temp array locations bracketing our temperature
    
    tdiff = abs(linetemps - layer%temp)
    
    
    Tlay1 = minloc(tdiff,1)
    
    
    if (linetemps(Tlay1) .lt. layer%temp) then
       Tlay2 = Tlay1+1
    else
       Tlay2 = Tlay1 - 1
    end if
    



    if (layer%temp .lt. linetemps(1)) then
       Tlay1 = 1
       Tlay2 = 2
    else if (layer%temp .gt. linetemps(nlinetemps)) then
       Tlay1 = nlinetemps - 1
       Tlay2 = nlinetemps
    endif


    
    ! get linear interpolation factor
    
    if (Tlay1 .gt. Tlay2) then
       torder = 1
       intfact = (log10(layer%temp) - log10(linetemps(Tlay2))) / (log10(linetemps(Tlay1)) - log10(linetemps(Tlay2)))
    else
       torder = 2
       intfact =  (log10(layer%temp) -log10(linetemps(Tlay1))) / (log10(linetemps(Tlay2)) - log10(linetemps(Tlay1)))
    endif

    if (layer%temp .gt. linetemps(nlinetemps)) then
       intfact =  (log10(layer%temp) -log10(linetemps(Tlay2))) / (log10(linetemps(Tlay2)) - log10(linetemps(Tlay1)))
    endif

    
    ! now get the files by their Tlay and pressure locations
    ! gases are identified by name (character!)
    ! this will allow the relevant linelists to be grabbed by T and P index
    
    
    
    ! loope through the gases
    
    totkappa = 0.0
    opd_lines = 0.0
    
    do i = 1, ngas
       if (layer%index .eq. 1) then
          write(*,*) "dropping in line opacities for ",trim(layer%gas(i)%name)
       end if
       
       write(lines1,"(A,A,A1,A,A1,I0,A1,I0)")"../LineLists/",trim(layer%gas(i)%name), &
            "/",trim(layer%gas(i)%name),"_",Tlay1,"_",layer%index
       write(lines2,"(A,A,A1,A,A1,I0,A1,I0)") "../LineLists/",trim(layer%gas(i)%name), &
            "/",trim(layer%gas(i)%name),"_",Tlay2,"_",layer%index
       

       
       
       
       ! now read in the line lists and sum, weighted by abudance/fraction
       
       iounit = index*10*i
       open(iounit,file=lines1, status='old')
!       write(*,*) "reading lines1 ", trim(lines1)
       read(iounit,*) kappa1
       close(iounit)
       do j = 1, nwave

          if (kappa1(j) .eq. 0.0) then
             logkap1(j) = -100.0
          else
             logkap1(j) = log10(kappa1(j))
          end if
       enddo
       open(iounit,file=lines2, status='old')
!       write(*,*) "reading lines2 ", trim(lines2)
       read(iounit,*) kappa2
       close(iounit)
       do j = 1,nwave
          if (kappa2(j) .eq. 0.0) then
             logkap2(j) = -100.0
          else
             logkap2(j) = log10(kappa2(j))
          end if
       enddo

       
       !  do the interpolation in log(cross section) !!!         !


       if (torder .eq. 1) then 
          logintkappa = ((logkap1 - logkap2)*intfact)+logkap2
       else if (torder .eq. 2) then
          logintkappa = ((logkap2 - logkap1)*intfact)+logkap1
       else
          write(*,*) "something wrong with interpolate order"           
       endif

       if (layer%temp .gt. linetemps(nlinetemps)) then
          logintkappa = ((logkap2 - logkap1)*intfact)+logkap2
       endif

       ! now we've got cross section for gas in layer - intkappa
       ! (cm2 / molecule)
       ! we want the optical depth.. 
       ! so we multiply by number density
       !(which is in /m3 to get column density ( / cm)
       ! and sum over the layer thickness (which we calculated in METRES)
       ! to get optical depth
       
       
       opd_lines = opd_lines + &
            ((layer%gas(i)%VMR * layer%ndens * layer%dz * 1d-4) &
            * (10.0**logintkappa))
    enddo
    
    
    
    
  end subroutine line_mixer
  
  
end module gas_mixing
