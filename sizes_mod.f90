  module sizes

  implicit none


  ! number of condensate species
  integer, parameter :: nclouds = 1

  ! number of gases for mixing
  integer, protected :: ngas

  ! number of patches
  integer,protected :: npatch
  
  ! declares and values size parameters for arrays 
  
  ! number of pressure layers
  integer, parameter :: nlayers = 64

  ! number of temperatures in line list grid
  integer, parameter :: nlinetemps = 27
 
  
  ! number of temperatures in CIA tables
  integer,parameter :: nciatemps = 198

  ! number of wavenumbers in lowres CIA tables
  integer,parameter :: ncwave = 1000

  
  ! number of radius bins for mie coeff files
  integer, parameter :: nrad = 40

  ! number of wavelengths in mie coeff files
  integer, parameter :: nmiewave = 180

  ! minimum radius for mie coeffs
  real, parameter :: rmin = 1e-5

  ! Volume ratio for sets in particle sizes
  real, parameter :: vrat = 2.2 

  ! max wave number
  integer, parameter :: maxwave = 19501
  
  ! number of wavelength/number bins in full range
  integer,protected :: nwave

  save

  
contains
  subroutine initwave(val)
    integer :: val
    nwave = val
  end subroutine initwave

  subroutine initgas(gasval)
    integer :: gasval
    ngas = gasval
  end subroutine initgas

  subroutine initpatch(pval)
    integer :: pval
    npatch =  pval
  end subroutine initpatch

  
end module sizes
