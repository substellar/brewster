module main
  
  implicit none

contains 
  subroutine forward(w1,w2,temp,logg,R2D2,gasname,molmass,logVMR,pcover,&
       do_clouds,cloudname,cloudrad,cloudsig,cloudprof,out_spec)
    
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
    
    
    ! input variables
    double precision,dimension(nlayers), INTENT(IN):: temp
    double precision,INTENT(IN) :: w1,w2
    real, INTENT(IN) :: R2D2, logg
    real,dimension(npatch) :: pcover
    integer,dimension(npatch):: do_clouds
    character(len=10),dimension(ngas),INTENT(IN) :: gasname
    double precision, dimension(ngas),INTENT(IN) :: molmass
    double precision, dimension(ngas,nlayers),INTENT(IN) :: logVMR
    character(len=10),dimension(npatch,nclouds), INTENT(IN) :: cloudname
    double precision, dimension(npatch,nlayers,nclouds),INTENT(IN) :: cloudrad
    double precision, dimension(npatch,nlayers,nclouds),INTENT(IN) :: cloudsig
    double precision, dimension(npatch,nlayers,nclouds),INTENT(IN) :: cloudprof
    double precision,dimension(2,nwave),INTENT(OUT) :: out_spec
    

    
    double precision, dimension(nlayers):: tmppress
    real:: metal,grav,test
    double precision,dimension(nwave) :: wdiff
    ! counters
    integer :: ch4index,ipatch,icloud,ilayer,iwave,igas,nw1,nw2
    real:: totcover, fboth, fratio, tstart,tfinish, opstart,opfinish
    real:: linstart, linfinish,distart, difinish

    
    ! HARD CODED to do line and CIA opacity calcs and everything
    ! apart from dust for first patch
    ! and copy to rest of patches before disort

    
    patch(1)%atm%temp = temp
    
    
      
    call set_pressure_scale

    ch4index=0
    do igas = 1, ngas
       
       
       patch(1)%atm%gas(igas)%name = gasname(igas)
       patch(1)%atm%gas(igas)%VMR = 10.**(logVMR(igas,:))
       patch(1)%atm%gas(igas)%molmass = molmass(igas)
       
       if (trim(patch(1)%atm(1)%gas(igas)%name) .eq. "ch4") then
          ch4index = igas
       endif

    end do
    
    
    grav = 10**(logg) / 100.
    
    
    ! now let's check how long this takes
    
    call cpu_time(tstart)
    
    
    
    ! now H2 and He fractions and mu for each layer
    
    do ilayer = 1, nlayers
       
       fboth = 1.0 - sum(patch(1)%atm(ilayer)%gas%VMR)
       
       ! hardcoded H/He ratio
       ! from solar abundance of 91.2 by number H, 8.7 by number He
       fratio = 0.84
       patch(1)%atm(ilayer)%fH2 = fratio * fboth
       patch(1)%atm(ilayer)%fHe = (1.0  - fratio) * fboth
       
       patch(1)%atm(ilayer)%mu = (patch(1)%atm(ilayer)%fH2 * XH2) + &
            (patch(1)%atm(ilayer)%fHe * XHe) + &
            sum(patch(1)%atm(ilayer)%gas%VMR * patch(1)%atm(ilayer)%gas%molmass)
       
    end do
    ! now we want the layer thickness in LENGTH units
    
    write(*,*) "Test line main L164. mu at layer 6 is: ", patch(1)%atm(6)%mu
    
    call layer_thickness(patch(1)%atm%press,patch(1)%atm%temp,grav,patch(1)%atm%dz)
    
    ! now get number density of layers
    ! number density in /m3:
    ! pressure is in mbar... so x100 to get N/m2
    
    patch(1)%atm%ndens = 100 * patch(1)%atm%press  / (K_BOLTZ * patch(1)%atm%temp)
    
    
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
    
    wdiff = abs(w1 - wavelen)
    nw1 = minloc(wdiff,1)
      
    wdiff = abs(w2 - wavelen)
    nw2 = minloc(wdiff,1)
    
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
    
    ! now put in the cloud details
    
    do ipatch = 1, npatch
       
       
       patch(ipatch)%cover = pcover(ipatch)
       
       patch(ipatch)%cloudy = do_clouds(ipatch)
       
       if (patch(ipatch)%cloudy) then
          do icloud = 1, nclouds
             patch(ipatch)%atm(1)%cloud(icloud)%name = cloudname(ipatch,icloud)
             patch(ipatch)%atm%cloud(icloud)%name = &
                  patch(ipatch)%atm(1)%cloud(icloud)%name
             patch(ipatch)%atm%cloud(icloud)%density = 10.**(cloudprof(ipatch,:,icloud))
             patch(ipatch)%atm%cloud(icloud)%rg = cloudrad(ipatch,:,icloud) 
             patch(ipatch)%atm%cloud(icloud)%rsig = cloudsig(ipatch,:,icloud)
             
          end do
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

    ! tests to see what's going to DISORT
    do ipatch = 1, npatch
       do ilayer = 1, nlayers
          !patch(ipatch)%atm(ilayer)%opd_lines = 0.
          !patch(ipatch)%atm(ilayer)%opd_CIA = 0.
       end do
    end do
    write(*,*) "Do clouds and cover: ", patch%cloudy, patch%cover, " OK?"

    ! now let put it all into DISORT
    
    write(*,*) "TEST: calling RUN_DISORT"
    write(*,*) " running between wavenum entries nw1, nw2: ", nw1,nw2
    
    call cpu_time(distart)
    
    call run_disort(out_spec(2,:),nw1,nw2)

    out_spec(1,:) = wavelen
    call cpu_time(difinish)
    
    
    call cpu_time(tfinish)
    
    write(*,*) "Time elapsed :", (tfinish - tstart), " seconds"
    
    write(*,*) "Opacity interpolations took : ", (opfinish - opstart), " seconds"
    
    write(*,*) "DISORT took : ", (difinish - distart), " seconds"
    
  end subroutine forward


end module main
