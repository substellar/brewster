module atmos_ops

  use sizes
  use phys_const
  
  implicit none


contains

  subroutine layer_thickness(press,temp,grav,dz)

    
    use sizes
    use phys_const
    
    implicit none
    
    !declare counters here to ensure they don't clash with counters from calling
    ! program
    integer:: i
    
    double precision :: p1,p2
    real, intent(IN) :: grav
    double precision, dimension(nlayers), intent(IN) :: press
    real, dimension(nlayers), intent(IN):: temp
    double precision, dimension(nlayers), intent(OUT):: dz
    
    
    ! The values in the pressure grid are geometric means for the layer. 
    ! Assume these are taken as sqrt(Pbot * Ptop)
    ! So (see derivation in notes)
    ! log(P1) = (3/2)*log(Pbar1) - (1/2)*log(Pbar2)
    ! log(P2) = 1/2(log (Pbar1) + log(Pbar2))
    
    ! we're using the hypsometric equation for this..
    
        
    do i = 1, nlayers-1 
       
       p1 = exp(((1.5)*log(press(i))) - ((0.5)*log(press(i+1))))
       p2 = exp((0.5)*(log(press(i) * press(i+1))))
       
       ! TK test line
      
       write(*,*) "TEST atmos_ops L45:  ", p1, p2

       ! no need to change P units to Pa, as we want the ratio
       
       dz(i)  = abs((R_GAS * temp(i) / grav) * log(p2 / p1))

        
       
    end do

    ! last layer needs special treatment
    ! is at the bottom of the atmosphere, so probably doesn't matter too much
    ! make it a mirror of layer above 

    dz(nlayers) = dz(nlayers - 1)
    
    
  end subroutine layer_thickness





  subroutine  set_temp_levels(layertemp,leveltemp)

    ! DISORT wants the temperature at levels, which we here retrieve from the
    ! Tbars that are in our atmospheric layers.
    ! This is kept in atmos_ops for consistency a
    integer :: i
    real, dimension(nlayers),INTENT(IN):: layertemp
    real,INTENT(OUT):: leveltemp(0:nlayers)


    leveltemp(0) = layertemp(1) - 0.5*(layertemp(2) - layertemp(1))

    ! if temperature is low and goes zero, lets just make the level equal to Tbar.
    
    if (leveltemp(0) .lt. 0) leveltemp(0) = layertemp(1)
    
    do i = 1, nlayers-1
       
       leveltemp(i) = layertemp(i) + 0.5*(layertemp(i+1) -layertemp(i))

       if (leveltemp(i) .lt. 0.0) leveltemp(i) = layertemp(i)
    end do


    leveltemp(nlayers) =  layertemp(nlayers) + 0.5*(layertemp(nlayers) -layertemp(nlayers-1))
    
  end subroutine set_temp_levels
  

  
  
end module atmos_ops
