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
  
  ! set up cloud and gas arrays, and atmosphere object
  type(a_gas) :: gas(ngas)
  type(a_cloud) :: cloud(ncloud)
  type(a_atmos) :: atm
  
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
  
!  if (tmppress .ne. press) then
!     write(*,*) "Input pressure scale doesn't match the line-list grid, please check and correct."
!     stop
!  endif

  ! set up the pressure layers - these are hard coded to match the line list
  ! pressure intervals to reduce processing

  atm%topdown = .true.
  atm%press=[1.e-02, 3.e-02, 1.e-1, 3.e-1,1.e0, 3.e0,1.e1, 3.e1, 1.e2, &
         3.e2,1.e3, 3.e3,1.e4, 3.e4,1.e5,3.e5] 
    
  atm%logP = log10(atm%press)


  
  ! Get VMRs, fixed in development for each gas, and write to all layers

  ! set gas order
  gas%topdown = .true.
  
  write(*,*) "Give VMR file"
  read(*,*) vmrfile

  open(10,file=vmrfile,status='old')

  
  do i = 1, ngas

     read(10,*) gas(i)%ID, fixVMR(i)

     gas(i)%VMR = fixVMR(i)
     ! TK test line
     write(*,*) gas(i)%VMR
  end do

  close(10)
  
  ! Also fudge the cloud params for development.
  ! Read in single values and copy to all layers

  !set cloud order
  cloud%topdown = .true.
  
  write(*,*) "Give cloud params file"
  read(*,*) cloudfile
  
  open(10,file=cloudfile,status="old")

  do i = 1, ncloud
     read(10,*) cloud(i)%id, tmp_n, tmp_rmax, tmp_rwidth
     cloud(i)%density = tmp_n(i)
     cloud(i)%rpeak = tmp_rmax(i)
     cloud(i)%rsigma = tmp_rwidth(i)
     ! TK test line
     write(*,*) cloud(i)%density
  end do


  write(*,*) "Give scale factor (R2/D2), log g(cgs)  and metallicity"
  read(*,*) R2D2, logg, metal

  ! put gravity in linear m/s units
  grav = 10**(logg) / 100.
  

  ! now we want the layer thickness in LENGTH units

  call layer_thickness(atm%press,atm%temp,grav,atm%dz)


  

  
  
end program main
