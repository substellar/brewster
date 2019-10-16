module common_arrays

  use sizes
  use define_types
  implicit none

  
  ! declares arrays that will be used alot
  
   ! the wavenumber array

   double precision, allocatable,dimension(:) :: wavenum, wavelen
   !real,dimension(nlinetemps) :: linetemps
   real,allocatable,dimension(:) :: press, linetemps
   integer, allocatable,dimension(:) :: gasnum
   integer,allocatable,dimension(:,:) :: cloudnum

   ! set up patchy atmosphere object, which is an array of patches
   type(a_patch),allocatable,dimension(:) :: patch

   save
 contains
   
   subroutine init_all
     if ( .NOT. allocated (wavenum)) allocate(wavenum(nwave))
     if ( .NOT. allocated (wavelen)) allocate(wavelen(nwave))
     if ( .NOT. allocated (gasnum)) allocate(gasnum(ngas))
     if ( .NOT. allocated (cloudnum)) allocate(cloudnum(npatch,nclouds))
     if ( .NOT. allocated (patch)) allocate(patch(npatch))     
     if ( .NOT. allocated (press)) allocate(press(nlayers))     
     if ( .NOT. allocated (linetemps)) allocate(linetemps(nlinetemps))     
     
   end subroutine init_all


     
end module common_arrays
