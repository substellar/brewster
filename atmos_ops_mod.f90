module atmos_ops

  use sizes
  use common_arrays
  use phys_const
  use define_types
  
  implicit none


contains

  subroutine set_pressure_scale


    integer:: i,ipatch

    do ipatch= 1, npatch
       patch(ipatch)%atm%index = (/ (i, i = 1, nlayers) /)
       patch(ipatch)%atm%press= press
       patch(ipatch)%atm%logP = log10(patch(ipatch)%atm%press)
    end do

  end subroutine set_pressure_scale


  
  subroutine layer_thickness(grav)

    
    implicit none
    
    !declare counters here to ensure they don't clash with counters from calling
    ! program
    integer:: i
    
    double precision :: p1, p2, R_spec
    real, intent(IN) :: grav
    
    ! The values in the pressure grid are geometric means for the layer. 
    ! Assume these are taken as sqrt(Pbot * Ptop)
    ! So (see derivation in notes)
    ! log(P1) = (3/2)*log(Pbar1) - (1/2)*log(Pbar2)
    ! log(P2) = 1/2(log (Pbar1) + log(Pbar2))
    
    ! we're using the hypsometric equation for this..
    ! careful! this requires SPECIFIC GAS CONSTANT
    
    
    do i = 1, nlayers
       
       R_spec = K_boltz / (patch(1)%atm(i)%mu * amu) 

       if (i .eq. nlayers) then
          p1 =  exp((0.5)*(log(patch(1)%atm(i-1)%press * patch(1)%atm(i)%press)))

          p2 = (patch(1)%atm(i)%press**2) / p1

          patch(1)%atm(i)%dp = p2 - p1

          ! no need to change P units to Pa, as we want the ratio
          
          patch(1)%atm(i)%dz  = abs((R_spec * patch(1)%atm(i)%temp / grav) * log(p1 /p2))


       else
          p1 = exp(((1.5)*log(patch(1)%atm(i)%press)) - &
               ((0.5)*log(patch(1)%atm(i+1)%press)))
          p2 = exp((0.5)*(log(patch(1)%atm(i)%press * patch(1)%atm(i+1)%press)))
          
          patch(1)%atm(i)%dp = p2 - p1
           
          ! no need to change P units to Pa, as we want the ratio
          
          patch(1)%atm(i)%dz  = abs((R_spec * patch(1)%atm(i)%temp / grav) * log(p1 /p2))
       end if
       
    end do
   

    
    ! last layer needs special treatment
    ! is at the bottom of the atmosphere, so probably doesn't matter too much
    ! make it a mirror of layer above 

    
  end subroutine layer_thickness





  subroutine  set_temp_levels(leveltemp)

    ! DISORT wants the temperature at levels, which we here retrieve from the
    ! Tbars that are in our atmospheric layers.
    ! This is kept in atmos_ops for consistency a
    integer :: i
    double precision :: p1,p2
    double precision,INTENT(OUT):: leveltemp(0:nlayers)

    ! we'll interpolate the level temps using Pbar, Tbar, and dP


    do i = 1, nlayers
   
       if (i .eq. 1) then
          
          p1 = exp(((1.5)*log(patch(1)%atm(i)%press)) - &
               ((0.5)*log(patch(1)%atm(i+1)%press)))
          p2 = exp((0.5)*(log(patch(1)%atm(i)%press * patch(1)%atm(i+1)%press)))

          leveltemp(i-1) = patch(1)%atm(i)%temp + &
               (((patch(1)%atm(i+1)%temp - patch(1)%atm(i)%temp) / &
               (patch(1)%atm(i+1)%logP  - patch(1)%atm(i)%logP)) * &
               (log10(p1) - patch(1)%atm(i)%logP))
          
          leveltemp(i) = patch(1)%atm(i)%temp + &
               (((patch(1)%atm(i+1)%temp - patch(1)%atm(i)%temp) / &
               (patch(1)%atm(i+1)%logP  - patch(1)%atm(i)%logP)) * &
               (log10(p2) - patch(1)%atm(i)%logP))
          
       elseif (i .eq. nlayers) then


          p1 =  exp((0.5)*(log(patch(1)%atm(i-1)%press * patch(1)%atm(i)%press)))          
          p2 = (patch(1)%atm(i)%press**2) / p1

          leveltemp(i) = patch(1)%atm(i)%temp + &
               (((patch(1)%atm(i)%temp - patch(1)%atm(i-1)%temp) / &
               (patch(1)%atm(i)%logP  - patch(1)%atm(i-1)%logP)) * &
               (log10(p2) - patch(1)%atm(i)%logP))

       else        

          p2 = exp((0.5)*(log(patch(1)%atm(i)%press * patch(1)%atm(i+1)%press)))

          leveltemp(i) = patch(1)%atm(i)%temp + &
               (((patch(1)%atm(i+1)%temp - patch(1)%atm(i)%temp) / &
               (patch(1)%atm(i+1)%logP  - patch(1)%atm(i)%logP)) * &
               (log10(p2) - patch(1)%atm(i)%logP))
       end if


       
    end do
    


    
  end subroutine set_temp_levels
  

  
  
end module atmos_ops
