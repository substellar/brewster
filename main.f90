program main

  use sizes
  use common_arrays
  use define_types
  use phys_const
  use atmos_ops
  use gas_mixing
  use cia_lowres
  use clouds
  use setup_disort
  
  implicit none


  ! some variables we need
  real, dimension(nlayers):: tmppress
  real:: metal,R2D2,logg, grav,test
  real,dimension(nwave) :: out_spec,patch_spec
  ! counters
  integer :: ch4index,ipatch,icloud,ilayer,iwave,igas
  real:: totcover, fboth, fratio, tstart,tfinish, opstart,opfinish
  real:: linstart, linfinish,distart, difinish
  
  ! some variables we'll get rid of later - just for development
  character(len=50) :: TPfile, VMRfile, cloudfile,fmt, junk
  character(len=10),dimension(ngas) :: gasname
  character(len=1):: response
  double precision, dimension(ngas) :: fixVMR
  double precision, dimension(nclouds) :: tmp_n, tmp_rmax, tmp_rwidth
  real:: isotemp
  logical:: do_clouds
  real :: w1, w2
  real, dimension(nwave) :: diff
  integer :: nw1, nw2
  

  ! set the wavelength range

  write(*,*) "give the wavelength range in microns"
  read(*,*) w1,  w2 
  
 

  
  ! get the temp profile and check that the pressure layers
  ! match the fixed values

  write(*,*) "Give the temperature profile"
  read(*,*) TPfile
  
  
  
  open(unit = 10, file=TPfile, status='old')

  do ilayer = 1, nlayers

     read(10,*) patch(1)%atm(ilayer)%temp

!     patch(:)%atm(i)%temp = patch(1)%atm(i)%temp
     
  end do

  close(10)
  

  call set_pressure_scale

  ! HARD CODED to do line and CIA opacity calcs and everything
  ! apart from dust for first patch
  ! and copy to rest of patches before disort


  ! Get VMRs, fixed in development for each gas, and write to all layers

  write(*,*) "Give VMR file"
  read(*,*) vmrfile

  fmt = "(A8,F14.12)"
  open(10,file=vmrfile,status='old')

  do igas = 1, ngas

     read(10,fmt) gasname(igas), fixVMR(igas)

     patch(1)%atm%gas(igas)%name = gasname(igas)
     patch(1)%atm%gas(igas)%VMR = fixVMR(igas)

      
     if (trim(patch(1)%atm(1)%gas(igas)%name) .eq. "ch4") then
        patch(1)%atm%gas(igas)%molmass = XCH4
        ch4index = igas
     endif
     if (trim(patch(1)%atm(1)%gas(igas)%name) .eq. "h2o") then
        patch(1)%atm%gas(igas)%molmass = XH2O
     end if
  end do
  

  close(10)
  

!  write(*,*) "Give scale factor (R2/D2), log g(cgs)  and metallicity"
!  read(*,*) R2D2, logg, metal

  ! put gravity in linear m/s units

  logg = 4.5
  grav = 10**(logg) / 100.


  ! now let's check how long this takes

  call cpu_time(tstart)


  
  ! now H2 and He fractions and mu for each layer

  do ilayer = 1, nlayers
     
     fboth = 1.0 - sum(patch(1)%atm(ilayer)%gas%VMR)

     ! hardcoded H/He ratio
     ! from solar abundance of 91.2 by number H, 8.7 by number He
     fratio = 0.84
     patch(1)%atm(ilayer)%fH2 = 0.84 * fboth
     patch(1)%atm(ilayer)%fHe = (1.0  - 0.84) * fboth

     patch(1)%atm(ilayer)%mu = (patch(1)%atm(ilayer)%fH2 * XH2) + &
          (patch(1)%atm(ilayer)%fHe * XHe) + &
          sum(patch(1)%atm(ilayer)%gas%VMR * patch(1)%atm(ilayer)%gas%molmass)
  
  end do
  ! now we want the layer thickness in LENGTH units

  write(*,*) "Test line main L164. mu at layer 6 is: ", patch(1)%atm(6)%mu
  
  call layer_thickness(patch(1)%atm%press,patch(1)%atm%temp,grav,patch(1)%atm%dz)

  write(*,*) "LAYER THICKNESS COMPLETE"

  ! now get number density of layers
  ! number density in /m3:
  ! pressure is in bar... so x1e5 to get N/m2
    
  patch(1)%atm%ndens = 1e5 * patch(1)%atm%press  / (K_BOLTZ * patch(1)%atm%temp)


  ! zero all the opacities
  do ipatch = 1, npatch
     do ilayer = 1, nlayers
        patch(ipatch)%atm(ilayer)%opd_scat = 0.0
        patch(ipatch)%atm(ilayer)%opd_CIA = 0.0
        patch(ipatch)%atm(ilayer)%opd_ext = 0.0
        patch(ipatch)%atm(ilayer)%opd_lines = 0.0
     end do
  end do

  ! get the wave array

  open(unit=10, file="../Linelists/wavegrid.dat", status='old')
  read(10,*) wavenum
  close(10)
  ! now we've got the wavenumber arrays we can set the index range for later

  wavelen = 1e4 / wavenum

  diff = abs(w1 - wavelen)
  nw1 = minloc(diff,1)
  
  diff = abs(w2 - wavelen)
  nw2 = minloc(diff,1)
 
  call cpu_time(opstart)
  call cpu_time(linstart)
  ! now mix the gases in each layer to get optical depth from lines
  do ilayer = 1, nlayers
     call line_mixer(patch(1)%atm(ilayer),patch(1)%atm(ilayer)%opd_lines,ilayer)
  end do

  call cpu_time(linfinish)

  write(*,*) "lines mixed in", (linfinish - linstart), "seconds. moving on to CIA"
  
  ! now let's get the CIA.  
  
  call get_cia_LR(grav,ch4index)

  
  call cpu_time(opfinish)
  
  ! COPY ALL PATCH 1 PROPERTIES TO OTHER PATCHES BEFORE DOING CLOUDS
  
  if (npatch .gt. 1) then
     do ipatch = 1, npatch
        patch(ipatch)%atm = patch(1)%atm
     end do
  end if

!  write(*,*) patch(1)%atm(5opd_lines(250000), patch(2)%atm(50)%opd_lines(250000)
!  write(*,*) patch(1)%atm(50)%opd_CIA(250000), patch(2)%atm(50)%opd_CIA(250000)

!  stop
   ! now put in the cloud details

  do ipatch = 1, npatch
     write(*,*) "Give covering fraction for patch", ipatch
     read(*,*) patch(ipatch)%cover

     write(*,*) "Is patch",ipatch," cloudy?"
     read(*,*) response
     
     patch(ipatch)%cloudy = .false.
     
     if(response .eq. 'y') then
        patch(ipatch)%cloudy = .true.
     endif
     
     if (patch(ipatch)%cloudy) then
        write(*,*) "Give cloud params file"
        read(*,*) cloudfile
        
        open(10,file=cloudfile,status="old")
        
        do icloud = 1, nclouds
           read(10,*) patch(ipatch)%atm(1)%cloud(icloud)%name
           do ilayer = 1, nlayers
              patch(ipatch)%atm(ilayer)%cloud(icloud)%name = &
                   patch(icloud)%atm(1)%cloud(icloud)%name
              read(10,*) patch(ipatch)%atm(ilayer)%cloud(icloud)%density,&
                   patch(ipatch)%atm(ilayer)%cloud(icloud)%rg, &
                   patch(ipatch)%atm(ilayer)%cloud(icloud)%rsig
           end do
        end do
        close(10)
        write(*,*) "calling cloud atlas"
        call cloudatlas(patch(ipatch)%atm)
     else
        do ilayer = 1, nlayers
           patch(ipatch)%atm(ilayer)%gg = 0.0
           patch(ipatch)%atm(ilayer)%opd_scat = 0.0
           patch(ipatch)%atm(ilayer)%opd_ext = 0.0
        end do
     end if

  end do
    
  totcover = sum(patch%cover)
  if (abs(1.0 - totcover) .gt. 0.001) then
     write(*,*) "Patch coverage doesn't sum to unity. Stopping."
     stop
  end if


  
  ! now let put it all into DISORT

  write(*,*) "TEST: calling RUN_DISORT"
  write(*,*) " running between wavenum entries nw1, nw2: ", nw1,nw2

  call cpu_time(distart)

  call run_disort(out_spec,nw1,nw2)

  


  call cpu_time(difinish)
  
  open(unit=20, file="test_spectrum.dat",status="new")

  write(20,*) "test spectrum in microns and W/m2/um(??)"
  
  do iwave = nw2, nw1
     write(20,*) wavelen(iwave), out_spec(iwave)
  end do
  
  close(20)
  
  call cpu_time(tfinish)
  
  write(*,*) "Time elapsed :", (tfinish - tstart), " seconds"
  
  write(*,*) "Opacity interpolations took : ", (opfinish - opstart), " seconds"
  
  write(*,*) "DISORT took : ", (difinish - distart), " seconds"
  
end program main
