module main
  
  implicit none

contains 
  subroutine forward(temp,logg,R2D2,gasname,ingasnum,molmass,logVMR,&
       pcover,do_clouds,incloudnum,cloudname,cloudrad,cloudsig,cloudprof,&
       inlinetemps,inpress,inwavenum,linelist,cia,ciatemp,use_disort,out_spec)
    
    use sizes
    use common_arrays
    use define_types
    use phys_const
    use atmos_ops
    use gas_opacity
    use clouds
    use setup_RT
    
    implicit none
    
    
    ! input variables
    double precision,dimension(nlayers), INTENT(IN):: temp
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
    integer,dimension(npatch,nclouds),intent(in):: incloudnum
    double precision,intent(inout) :: inwavenum(:)
    real,dimension(nlinetemps) :: inlinetemps
    real,dimension(nlayers) :: inpress
    double precision,intent(inout):: linelist(:,:,:,:)
    real, dimension(nciatemps), intent(in):: ciatemp
    real, intent(inout) :: cia(:,:,:)
    integer, dimension(ngas), intent(in) :: ingasnum
    double precision,allocatable, dimension(:,:),INTENT(OUT) :: out_spec
    double precision,allocatable, dimension(:)::specflux
    real:: metal,grav,test
    double precision,allocatable :: wdiff(:)
    ! counters
    integer :: ch4index,ipatch,icloud,ilayer,iwave,igas,nw1,nw2,use_disort
    real:: totcover, fboth, fratio, tstart,tfinish, opstart,opfinish
    real:: linstart, linfinish,distart, difinish
    logical :: disorting

    ! Are we using DISORT
    disorting = use_disort
    ! HARD CODED to do line and CIA opacity calcs and everything
    ! apart from dust for first patch
    ! and copy to rest of patches before disort
    allocate(wdiff(nwave))
    call init_wscales
    !wavenum[nwave],wavelen[nwave])


    do ipatch = 1, npatch
       call init_column(patch(ipatch)%atm)
    end do
    
    patch(1)%atm%temp = temp
   
    wavenum = inwavenum
    linetemps = inlinetemps
    press = inpress
    gasnum = ingasnum
    cloudnum = incloudnum

    call set_pressure_scale

    ! TEST LINE - set artificial temp profile
    !patch(1)%atm%temp = 300. + (((patch(1)%atm%press/1000.) / 0.1)**1.1)


    
    ch4index=0
    do igas = 1, ngas
       
       
       patch(1)%atm%gas(igas)%name = gasname(igas)
       patch(1)%atm%gas(igas)%num = gasnum(igas)       
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
    
!    write(*,*) "Test line main L164. mu at layer 6 is: ", patch(1)%atm(6)%mu
    call layer_thickness(grav)
    
    ! now get number density of layers
    ! number density in /m3:
    ! pressure is in bar... so x1e5 to get N/m2
    
    patch(1)%atm%ndens = 1.e+5 * patch(1)%atm%press  / (K_BOLTZ * patch(1)%atm%temp)
    
      ! zero all the opacities
    do ipatch = 1, npatch
       do ilayer = 1, nlayers
          patch(ipatch)%atm(ilayer)%opd_scat = 0.0
          patch(ipatch)%atm(ilayer)%opd_CIA = 0.0
          patch(ipatch)%atm(ilayer)%opd_ext = 0.0
          patch(ipatch)%atm(ilayer)%opd_lines = 0.0
          patch(ipatch)%atm(ilayer)%opd_rayl = 0.0
          
       end do
    end do
    
    wavelen = 1e4 / wavenum
    
 !   wdiff = abs(w1 - wavelen)
 !   nw1 = minloc(wdiff,1)
      
 !   wdiff = abs(w2 - wavelen)
 !   nw2 = minloc(wdiff,1)
    
    call cpu_time(opstart)
    ! now mix the gases in each layer to get optical depth from lines
    do ilayer = 1, nlayers
       call line_mixer(patch(1)%atm(ilayer),patch(1)%atm(ilayer)%opd_lines,ilayer,linelist)
 
       ! now the Rayleigh scattering
       call get_ray(ilayer,ch4index)
       
    end do

    ! now let's get the CIA.  
    call get_cia(cia,ciatemp,grav,ch4index)

    
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
             ! cloud profile given in units of log(n_cond/n_gas)
             ! now convert to density
             patch(ipatch)%atm%cloud(icloud)%density = &
                  patch(1)%atm%ndens * (10.**cloudprof(ipatch,:,icloud))
             patch(ipatch)%atm%cloud(icloud)%rg = cloudrad(ipatch,:,icloud) 
             patch(ipatch)%atm%cloud(icloud)%rsig = cloudsig(ipatch,:,icloud)
             
          end do ! cloud loop

          
          ! in case of simple/generic/mixed cloud we won't be doing Mie coeffs
          ! we'll just use the density, rg, and rsig as dtau, w0 and gg
          ! for the cloud

          !if (trim(patch(ipatch)%atm(1)%cloud(1)%name) .eq. 'mixto') then
          if (cloudnum(ipatch,1) .eq. 99) then
             if (nclouds .ne. 1) then
                write(*,*) "Error: mixto cloud case should have nclouds = 1"
                stop
             else
                ! FOR GREY SIMPLE CASE (MIXTO): put DTAU_cloud in cloudprofile
                ! and albedo in rg, and asymmetry in rsig
                ! There must only be one cloud for this case
                do ilayer= 1, nlayers
                   patch(ipatch)%atm(ilayer)%opd_ext = &
                        cloudprof(ipatch,ilayer,1)
                   patch(ipatch)%atm(ilayer)%opd_scat = &
                        patch(ipatch)%atm(ilayer)%opd_ext * &
                        patch(ipatch)%atm(ilayer)%cloud(1)%rg
                   patch(ipatch)%atm(ilayer)%gg = &
                        patch(ipatch)%atm(ilayer)%cloud(1)%rsig
                end do ! layer loop
             end if
          else 
             write(*,*) "calling cloud atlas"
             call cloudatlas(patch(ipatch)%atm)
          end if
       else
          do ilayer = 1, nlayers
             patch(ipatch)%atm(ilayer)%gg = 0.0
             patch(ipatch)%atm(ilayer)%opd_scat = 0.0
             patch(ipatch)%atm(ilayer)%opd_ext = 0.0
          end do
       end if
       
    end do ! patch loop
    
    totcover = sum(patch%cover)
    if (abs(1.0 - totcover) .gt. 0.001) then
       write(*,*) "Patch coverage doesn't sum to unity. Stopping."
       stop
    end if

    ! tests to see what's going to DISORT
    !do ipatch = 1, npatch
    !   do ilayer = 1, nlayers
    !      patch(ipatch)%atm(ilayer)%opd_lines = 0.
    !      patch(ipatch)%atm(ilayer)%opd_CIA = 0.
    !       end do
    !    end do

    ! test line
    call cpu_time(distart)

    
    call run_RT(specflux,disorting)

    allocate(out_spec(2,nwave))
    out_spec(1,:) = wavelen
    out_spec(2,:) = specflux * R2D2
    ! scale by r2d2
    
    call cpu_time(difinish)
    
    
    call cpu_time(tfinish)
    
    write(*,*) "Time elapsed :", (tfinish - tstart), " seconds"
    
    write(*,*) "Opacity interpolations took : ", (opfinish - opstart), " seconds"
    
    write(*,*) "RT took : ", (difinish - distart), " seconds"

    deallocate(wavelen,wavenum)
    do ipatch = 1, npatch
       do ilayer= 1, nlayers
          deallocate(patch(ipatch)%atm(ilayer)%opd_ext,&
               patch(ipatch)%atm(ilayer)%opd_scat,&
               patch(ipatch)%atm(ilayer)%opd_lines, &
               patch(ipatch)%atm(ilayer)%opd_cia, &
               patch(ipatch)%atm(ilayer)%gg)
       end do
    end do
   
  end subroutine forward

  
end module main
