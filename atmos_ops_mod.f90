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
    
    real :: p1,p2
    real, intent(IN) :: grav
    real, dimension(nlayers), intent(IN) :: press,temp
    real, dimension(nlayers), intent(OUT):: dz
    
    
    ! The values in the pressure grid are geometric means for the layer. 
    ! Assume these are taken as sqrt(Pbot * Ptop)
    ! So (see derivation in notes)
    ! log(P1) = (3/2)*log(Pbar1) - log(Pbar2)
    ! log(P2) = 1/2(log (Pbar1) + log(Pbar2))
    
    ! we're using the hypsometric equation for this..I
    
        
    do i = 1, nlayers-1 
       
       p1 = exp(((3/2)*log(press(i))) - log(press(i+1)))
       p2 = exp((1/2)*(log(press(i)) + log(press(i+1))))
       
       ! TK test line
       write(*,*) p1, p2
       
       
       dz(i)  = abs((R_GAS * temp(i) / grav) * log(p2 / p1))
       
    end do

    ! last layer needs special treatment
    ! is at the bottom of the atmosphere, so probably doesn't matter too much
    ! make it a mirror of layer above 

    dz(nlayers) = dz(nlayers - 1)
    
    
  end subroutine layer_thickness
  
  
  
  
  
end module atmos_ops
