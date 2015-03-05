module phys_const

  implicit none

  double precision, parameter :: pi = 3.14159274d0

  ! light in m/s
  real, parameter :: c = 299792458.00

  !   unknown to the Romans 
  double precision,  parameter :: ZERO = 0.d0 

  !   universal gas constant (J/mol/K)
  double precision, parameter :: R_GAS = 8.3144621 
  
  !   Stefan-Boltzmann constant (W/m^2/K^4)
  double precision, parameter :: STEFBOLTZ = 5.670373d-8

  !   Avogadros number (#/mol)
  double precision, parameter :: AVOGADRO = 6.02d23 

  !   Boltzmann constant (J/K)
  double precision, parameter :: K_BOLTZ = R_GAS / AVOGADRO 


  ! Atomic mass unit (kg)
  double precision, parameter :: amu = 1.6605402d-27
  
  
  save




end module phys_const
