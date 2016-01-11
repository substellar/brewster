module common_arrays

  use sizes
  use define_types
  implicit none

  
  ! declares arrays that will be used alot
  
   ! the wavenumber array

   double precision, allocatable,dimension(:) :: wavenum, wavelen
   real,dimension(nlinetemps) :: linetemps
   real,dimension(nlayers) :: press
   integer, dimension(ngas) :: gasnum

   save
 contains
   subroutine init_wscales
     if ( .NOT. allocated (wavenum)) allocate(wavenum(nwave))
     if ( .NOT. allocated (wavelen)) allocate(wavelen(nwave))
     
   end subroutine init_wscales
     
     
end module common_arrays
