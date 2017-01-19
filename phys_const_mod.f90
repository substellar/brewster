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

  ! Boltzmann constant (cgs)
  double precision, parameter :: Kbolt_cgs = 1.38064853d-16

  ! Atomic mass unit (kg)
  double precision, parameter :: amu = 1.6605402d-27
  
  ! some molecular masses
   double precision, parameter ::   XH2 = 2.01588d0
   double precision, parameter ::   XHe=4.002602d0
   double precision, parameter ::   XH2O=18.01528d0
   double precision, parameter ::   XCH4=16.04276d0
   double precision, parameter ::   XCO=28.0104d0
   double precision, parameter ::   XNH3=17.03056d0
   double precision, parameter ::   XPH3=33.997582d0
   double precision, parameter ::   XH2S=65.13994d0
   double precision, parameter ::   XTiO=63.8794d0
   double precision, parameter ::   XVO=66.9409d0
   double precision, parameter ::   XFeH=56.85494d0
   double precision, parameter ::   XFe=55.84d0
   double precision, parameter ::   XCrH=53.00404d0
   double precision, parameter ::   XH=1.008d0
   double precision, parameter ::   XN2=28.0134d0
   double precision, parameter ::   XCO2=44.0d0


  
  save




end module phys_const
