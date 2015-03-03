module phys_const_mod

  implicit none

  double precision, parameter :: pi = 3.14159274d0

  ! light in cm/s
  real, parameter :: c = 29979245800.00

  !   unknown to the Romans 
  double precision,  parameter :: ZERO = 0.d0 

  !   universal gas constant (erg/mol/K)
  double precision, parameter :: R_GAS = 8.3143d7 
  
  !   Stefan-Boltzmann constant (erg/cm^2/s/K^4)
  double precision, parameter :: STEFBOLTZ = 5.67d-5

  !   Avogadros number (#/mol)
  double precision, parameter :: AVOGADRO = 6.02d23 

  !   Boltzmann constant (erg/K)
  double precision, parameter :: K_BOLTZ = R_GAS / AVOGADRO 


  ! Atomic mass unit (grams)
  double precision, parameter :: amu = 1.6605402d-24
  
  
  save




end module phys_const_mod
