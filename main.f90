program main

  implicit none

  use sizes_mod
  use common_arrays_mod



  ! some variables
  real, dimension(nlayers):: tmppress



  ! some variables we'll get rid of later - just for development
  character(len=50) :: TPfile, VMRfile, dustfile
  double precision, dimension(ngas) :: fixVMR
  
  ! get the temp profile and check that the pressure layers match the fixed values

  
  write(*,*) "give the T-P profile file"
  read(*,*) TPfile

  open(unit = 10, file=TPfile, status='old')

  do i = 1, nlayer

     read(10,*) tmppress(i), temp(i)

  end do

  close(10)
  
  if (tmppress .ne. press) then
     write(*,*) "Input pressure scale doesn't match the line-list grid, please check and correct."
     stop
  endif

  ! Get VMRs, fixed in development for each gas, and write to all layers
  
  write(*,*) "Give VMR file"
  read(*,*) vmrfile

  open(10,file=vmrfile,status='old')

  
  do i = 1, ngas

     read(10,*) gasid(i), fixVMR(i)

     VMR(i,:) = fixVMR(i)
  end do

  close(10)
  

  

  
end program main
