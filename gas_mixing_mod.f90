module gas_mixing



  implicit none


contains

  subroutine line_mixer(layer,opd_lines,index)

    use sizes
    use phys_const
    use define_types
    use common_arrays


    implicit none
    
    type(a_layer), intent(IN) :: layer
    double precision, dimension(nwave), intent(OUT) :: opd_lines
    integer :: Tlay1, Tlay2, torder,index
    double precision, dimension(nwave) :: kappa1,kappa2,logintkappa,totkappa
    double precision, dimension(nwave) :: logkap1, logkap2
    real, dimension(nlinetemps):: tdiff
    double precision :: ndens, intfact, junk
    character(len=50) :: lines1,lines2,name

    ! counters
    integer:: i, j
    
    ! set up the line temp array - check this is working properly
    
    call set_line_temps
    
    ! TK test line
    !write(*,*) "Line temperatures: ", linetemps
    
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
       ! TK test line
       write(*,*) "torder = ",torder
       write(*,*) linetemps(Tlay2),linetemps(Tlay1)
    else
       torder = 2
       intfact =  (log10(layer%temp) -log10(linetemps(Tlay1))) / (log10(linetemps(Tlay2)) - log10(linetemps(Tlay1)))
       ! TK test line
       write(*,*) "torder = ",torder
       write(*,*) linetemps(Tlay2),linetemps(Tlay1)
    endif

    if (layer%temp .gt. linetemps(nlinetemps)) then
       intfact =  (log10(layer%temp) -log10(linetemps(Tlay2))) / (log10(linetemps(Tlay2)) - log10(linetemps(Tlay1)))
    endif

    
    ! now get the files by their Tlay and pressure locations
    ! gases are identified by name (character!)
    ! this will allow the relevant linelists to be grabbed by T and P index
    
    
    
    ! loope through the gases
    
    totkappa = 0.0
    
    do i = 1, ngas
       
       write(lines1,"(A,A,A1,A,A1,I0,A1,I0)")"../LineLists/",trim(layer%gas(i)%name), &
            "/",trim(layer%gas(i)%name),"_",Tlay1,"_",layer%index
       write(lines2,"(A,A,A1,A,A1,I0,A1,I0)") "../LineLists/",trim(layer%gas(i)%name), &
            "/",trim(layer%gas(i)%name),"_",Tlay2,"_",layer%index
       
       ! TK test line
       write(*,*) lines1
       write(*,*) lines2
       
       
       
       ! now read in the line lists and sum, weighted by abudance/fraction
       
       
       open(15,file=lines1, status='old')
       do j =1, listheadlines
          read(15,*)
       enddo
       do j = 1, nwave
          read(15,*) wavenum(j), kappa1(j)
          if (kappa1(j) .eq. 0.0) then
             logkap1(j) = -100.0
          else
             logkap1(j) = log10(kappa1(j))
          end if
       enddo
       close(15)
       open(15,file=lines2, status='old')
       do j =1, listheadlines
          read(15,*)
       enddo
       do j = 1,nwave
          read(15,*) junk, kappa2(j)
          if (kappa2(j) .eq. 0.0) then
             logkap2(j) = -100.0
          else
             logkap2(j) = log10(kappa2(j))
          end if
       enddo
       close(15)
       
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
       
       totkappa = totkappa + (layer%gas(i)%VMR * (10**logintkappa))
    enddo
    
    ! now we've got total cross section for layer - totkappa (cm2 / molecule)
    ! we want the optical depth.. 
    ! so we multiply by number density to get column density ( / cm)
    ! and sum over the layer thickness (which we calculated in METRES
    
    
    ! optical depth
    
    ! something here is making opd_line zero!
    opd_lines = totkappa * layer%ndens * layer%dz  * 10**(-4.0) 
    
    
    ! TK test output
    write(*,*) "line mixer L142 test.. layer%temp =", layer%temp
    write(*,*) "line mixer L142 test.. layer%dz =", layer%dz
    write(*,*) "line mixer L142 test.. ndens =", layer%ndens
    write(*,*) "line mixer L142 test.. totkappa 1000 =", totkappa(1000)
    
    
    
    
  end subroutine line_mixer
  
  
end module gas_mixing
