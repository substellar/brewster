module common_arrays

  use sizes
  use define_types
  implicit none

  
  ! declares arrays that will be used alot
  
   ! the wavenumber array

   double precision, dimension(nwave) :: wavenum, wavelen
   real,dimension(nlinetemps) :: linetemps
   real,dimension(nlayers) :: press
   integer, dimension(ngas) :: gasnum

   save


end module common_arrays
