module sizes

  implicit none

  ! declares and values size parameters for arrays 
  
  ! number of pressure layers
  integer, parameter :: nlayers = 16

  ! number of wavelength/number bins in full range
  integer, parameter :: nwave = 319188


  ! number of gases for mixing
  integer, parameter :: ngas = 5


  ! number of condensate species
  integer, parameter :: ncloud = 1

  
  save

  
end module sizes
