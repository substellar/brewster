program main

  use sizes
  use common_arrays
  use define_types
  use phys_const
  use atmos_ops
  
  implicit none


  ! some variables we need
  real, dimension(nlayers):: tmppress
  real:: metal,R2D2,logg, grav,test
  
  ! counters
  integer :: i
  
  
  ! some variables we'll get rid of later - just for development
  character(len=50) :: TPfile, VMRfile, cloudfile
  double precision, dimension(ngas) :: fixVMR
  double precision, dimension(ncloud) :: tmp_n, tmp_rmax, tmp_rwidth



  ! get the temp profile and check that the pressure layers
  ! match the fixed values

  
  write(*,*) "give the T-P profile file"
!  read(*,*) TPfile

!  open(unit = 10, file=TPfile, status='old')

!  do i = 1, nlayer

!     read(10,*) tmppress(i), atm%temp(i)

  
!  end do

!  close(10)
  
!  if (tmppress .ne. atm%press) then
!     write(*,*) "Input pressure scale doesn't match the line-list grid, please check and correct."
!     stop
!  endif

  call set_pressure_scale

  ! TK test line
  write(*,*) atm%press
  
  ! Get VMRs, fixed in development for each gas, and write to all layers

  ! set gas order
!  gas%topdown = .true.
  
  write(*,*) "Give VMR file"
  read(*,*) vmrfile

  open(10,file=vmrfile,status='old')

  
  do i = 1, ngas

     read(10,*) atm%gas(i)%name, fixVMR(i)

     atm%gas(i)%VMR = fixVMR(i)
     ! TK test line
     write(*,*) atm%gas(i)%VMR
  end do

  close(10)
  
  ! Also fudge the cloud params for development.
  ! Read in single values and copy to all layers

  !set cloud order
!  cloud%topdown = .true.
  
  write(*,*) "Give cloud params file"
  read(*,*) cloudfile
  
  open(10,file=cloudfile,status="old")

  do i = 1, ncloud
     read(10,*) atm%cloud(i)%name, tmp_n, tmp_rmax, tmp_rwidth
     atm%cloud(i)%density = tmp_n(i)
     atm%cloud(i)%rpeak = tmp_rmax(i)
     atm%cloud(i)%rsigma = tmp_rwidth(i)
     ! TK test line
     write(*,*) atm%cloud(i)%density
  end do


  write(*,*) "Give scale factor (R2/D2), log g(cgs)  and metallicity"
  read(*,*) R2D2, logg, metal

  ! put gravity in linear m/s units
  grav = 10**(logg) / 100.
  

  ! now we want the layer thickness in LENGTH units

  call layer_thickness(atm%press,atm%temp,grav,atm%dz)


  ! now mix the gases in each layer to get optical depth from lines


  
  

  

  
  
end program main
