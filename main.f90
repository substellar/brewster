program main

  use sizes
  use common_arrays
  use define_types
  use phys_const
  use atmos_ops
  use gas_mixing
  use setup_disort
  
  implicit none


  ! some variables we need
  real, dimension(nlayers):: tmppress
  real:: metal,R2D2,logg, grav,test
  real,dimension(nwave) :: dis_spec
  ! counters
  integer :: i
  
  
  ! some variables we'll get rid of later - just for development
  character(len=50) :: TPfile, VMRfile, cloudfile,fmt, junk
  character(len=10),dimension(ngas) :: gasname
  character(len=1):: response
  double precision, dimension(ngas) :: fixVMR
  double precision, dimension(ncloud) :: tmp_n, tmp_rmax, tmp_rwidth
  real:: isotemp
  logical:: do_clouds
  real :: w1, w2
  real, dimension(nwave) :: diff
  integer :: nw1, nw2
    
  

  ! set the wavelength range
  w1 = 1.25
  w2 = 1.75
  
 

  
  ! get the temp profile and check that the pressure layers
  ! match the fixed values



  write(*,*) "Are we doing clouds?"
  read(*,*) response

  do_clouds = .false.
  
  if(response .eq. 'y') then
     do_clouds = .true.
  endif

  
  write(*,*) "Give the temperature profile"
  read(*,*) TPfile


  atm%temp = isotemp

  open(unit = 10, file=TPfile, status='old')

 do i = 1, nlayers

     read(10,*) atm(i)%temp

  
  end do

  close(10)
  
!  if (tmppress .ne. atm%press) then
!     write(*,*) "Input pressure scale doesn't match the line-list grid, please check and correct."
!     stop
!  endif

  call set_pressure_scale

  ! TK test line
  write(*,*) "MAIN L71 TEST:", atm%press
  
  ! Get VMRs, fixed in development for each gas, and write to all layers

  ! set gas order
!  gas%topdown = .true.
  
  write(*,*) "Give VMR file"
  read(*,*) vmrfile

  fmt = "(A8,F6.4)"
  open(10,file=vmrfile,status='old')

  do i = 1, ngas

     read(10,fmt) gasname(i), fixVMR(i)

!     write(gasname(i),*) junk(1:8)
!     read(junk(9:14),*) fixVMR(i)
     atm%gas(i)%name = gasname(i)
     atm%gas(i)%VMR = fixVMR(i)
  end do

  close(10)
  
  ! Also fudge the cloud params for development.
  ! Read in single values and copy to all layers

  !set cloud order
!  cloud%topdown = .true.
  
  if (do_clouds) then
     write(*,*) "Give cloud params file"
     read(*,*) cloudfile
     
     open(10,file=cloudfile,status="old")
  
     do i = 1, ncloud
        read(10,*) atm%cloud(i)%name, tmp_n, tmp_rmax, tmp_rwidth
        atm%cloud(i)%density = tmp_n(i)
        atm%cloud(i)%rpeak = tmp_rmax(i)
        atm%cloud(i)%rsigma = tmp_rwidth(i)
        ! TK test line
        write(*,*) "MAIN L116 TEST: ", atm%cloud(i)%density
     end do
  end if

!  write(*,*) "Give scale factor (R2/D2), log g(cgs)  and metallicity"
!  read(*,*) R2D2, logg, metal

  ! put gravity in linear m/s units

  !grav = 10**(logg) / 100.
  grav = 10.

  ! now we want the layer thickness in LENGTH units
  ! TK TEST
  write(*,*) "TEST: Calling layer_thickness"
  
  call layer_thickness(atm%press,atm%temp,grav,atm%dz)

  write(*,*) "LAYER THICKNESS COMPLETE"

  ! now get number density of layers
      ! number density in /m3:
    
  atm%ndens = atm%press  / (K_BOLTZ * atm%temp)




  ! now mix the gases in each layer to get optical depth from lines
  do i = 1, nlayers
     call line_mixer(atm(i),atm(i)%opd_lines,i)
  end do

  
  ! now we've got the wavenumber arrays we can set the index range for later

  wavelen = 1e4 / wavenum

  diff = abs(w1 - wavelen)
  nw1 = minloc(diff,1)
  
  diff = abs(w2 - wavelen)
  nw2 = minloc(diff,1)
 

  
  
  ! zero all the missing bits
  do i = 1, nlayers
     atm(i)%opd_scat = 0.0
     atm(i)%opd_CIA = 0.0
  end do

  ! add up all the taus to get extinction
  do i = 1, nlayers
     atm(i)%opd_ext = atm(i)%opd_scat + atm(i)%opd_lines + atm(i)%opd_CIA
  end do


  write(junk,"(A,I0,A)") "test_line_opacities_layer_",8,".txt"
  open(unit=20,file=junk,status="new")
  write(20,*) "written at line 180 main"
  do i = 1, nwave
     write(20,*) wavelen(i), atm(8)%opd_ext(i)
  end do
  close(20)


  ! scattering opd stuff needed here and CIA!!!


  
  ! now let put it all into DISORT

  write(*,*) "TEST: calling RUN_DISORT"
  write(*,*) " running between wavenum entries nw1, nw2: ", nw1,nw2

  call run_disort(dis_spec,do_clouds,nw1,nw2)

  ! TK test edit Just a subset of wavelength space... say 1 - 5um

  

 
  open(unit=20, file="test_spectrum.dat",status="new")
  write(20,*) "test spectrum in microns and W/m2/um(??)"
  do i = nw2, nw1
     write(20,*) wavelen(i), dis_spec(i)
  end do
     
  close(20)
  
  
end program main
