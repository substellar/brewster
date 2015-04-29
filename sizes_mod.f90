module sizes

  implicit none

  ! declares and values size parameters for arrays 
  
  ! number of pressure layers
  integer, parameter :: nlayers = 61

  ! number of temperatures in line list grid
  integer, parameter :: nlinetemps = 27
  
  ! number of wavelength/number bins in full range
  integer, parameter :: nwave = 319188


  ! number of gases for mixing
  integer, parameter :: ngas = 2


  ! number of condensate species
  integer, parameter :: ncloud = 1

  ! number of lines of header in linelists
  integer,parameter :: listheadlines = 24

  ! number of temperatures in CIA tables
  integer,parameter :: nciatemps = 119
  
  save

  
end module sizes
