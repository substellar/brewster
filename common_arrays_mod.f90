module common_arrays

  use sizes
  use define_types
  implicit none

  
  ! declares arrays that will be used alot
  
  ! set up cloud and gas arrays, and atmosphere object
   type(a_layer) :: atm(nlayers)

  save

  ! set up the pressure layers - these are hard coded to match the line list
  ! pressure intervals to reduce processing
contains
  subroutine set_pressure_scale

    integer:: i
    
    atm%index = (/ (i, i = 1, nlayers) /)

    atm%press=[1.e-02, 3.e-02, 1.e-1, 3.e-1,1.e0, 3.e0,1.e1, 3.e1, 1.e2, &
         3.e2,1.e3, 3.e3,1.e4, 3.e4,1.e5,3.e5] 
    
    atm%logP = log10(atm%press)


  end subroutine set_pressure_scale
  
  subroutine set_line_temps

    real,dimension(nlinetemps) :: linetemps

    linetemps = [75., 85., 100., 120., 140., 160., 180., 200., 230., 260., &
           300., 350., 400.,500., 650., 800., 1000., 1200., 1400., 1600.,&
           1800., 2000., 2300., 2600.,3000., 3500., 4000.]


  end subroutine set_line_temps



end module common_arrays
