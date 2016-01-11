module gas_mixing



contains

  subroutine line_mixer(layer,opd_lines,index,linelist)

    use sizes
    use common_arrays
    use phys_const
    use define_types


    implicit none

    
    type(a_layer),intent(inout) :: layer
    double precision, intent(inout) :: opd_lines(:)
    double precision,intent(in):: linelist(:,:,:,:)
    integer :: Tlay1, Tlay2, torder, index, iounit
    double precision, allocatable,dimension(:,:) :: logkap1, logkap2,logintkappa
    real, dimension(nlinetemps):: tdiff
    double precision :: ndens, intfact, junk
    character(len=50) :: lines1,lines2,name
    
    ! counters
    integer:: igas, j
    
    allocate(logkap1(ngas,nwave),logkap2(ngas,nwave),logintkappa(ngas,nwave))
    
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

    
    
    opd_lines = 0.0
    

    !! CODE FUCKS UP HERE.. segfault 
    logkap1 = linelist(:,layer%index,Tlay1,:)
    logkap2 = linelist(:,layer%index,Tlay2,:) 
       
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

       
    do igas = 1, ngas       
       opd_lines = opd_lines + &
            ((layer%gas(igas)%VMR * layer%ndens * layer%dz * 1d-4) &
            * (10.0**logintkappa(igas,:)))
    enddo
    
    
    deallocate(logkap1, logkap2,logintkappa)
    
  end subroutine line_mixer
  
  
end module gas_mixing
