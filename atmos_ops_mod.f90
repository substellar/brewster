module atmos_ops

  use sizes
  use phys_const
  
  implicit none


contains

  subroutine layer_thickness(press,temp,gravity,dz)

    
    use sizes
    use phys_const
    
    implicit none
    
    !declare counters here to ensure they don't clash with counters from calling
    ! program
    integer:: i
    
    real :: p1,p2
    real, intent(IN) :: gravity
    real, dimension(nlayers), intent(IN) :: press,temp
    real, dimension(nlayers), intent(OUT):: dz
    
    
    ! We are assuming that the values in press(nlayers) are the pressures at
    ! the middle of each layer, so need to calculate the top and bottom pressures
    ! first
    
    ! first and last layers need special treatment
    ! lets make them symetrical about central pressure
    
    p1 = ((press(2) - press(1)) / 2) + press(1)
    p2 = ((press(2) - press(1)) / 2) - press(1)
    
    ! we're using the hypsometric equation for this..Is this OK?!!
    
    dz(1)  = abs((R_GAS * temp(1) / gravity) * log(p2 / p1))
    
    p1 = ((press(nlayers) - press(nlayers-1)) / 2) + press(1)
    p2 = ((press(nlayers) - press(nlayers-1)) / 2) - press(1)
    
    dz(nlayers)  = abs((R_GAS * temp(nlayers) / gravity) * log(p2 / p1))
    
    
    
    do i = 2, nlayers-1 
       
       p1 = ((press(i-1) - press(i)) / 2) + press(i)
       p2 =((press(i+1) - press(i)) / 2) + press(i)
       
       ! TK test line
       write(*,*) p1, p2
       
       
       dz(i)  = abs((R_GAS * temp(i) / gravity) * log(p2 / p1))
       
    end do
    
    
    
  end subroutine layer_thickness
  
  
  
  
  
end module atmos_ops
