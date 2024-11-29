module main
  
  implicit none

contains 
  subroutine forward(temp,logg,R2D2,gasname,molmass,logVMR,&
       pcover,do_clouds,incloudnum,cloudname,cloudrad,cloudsig,cloudprof,&
       inlinetemps,inpress,inwavenum,linelist,cia,ciatemp,use_disort,clphot,&
       othphot,do_cf,do_bff,bff,out_spec,clphotspec,othphotspec,cf)
    
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
    double precision,INTENT(IN):: temp(:)
    real, INTENT(IN) :: R2D2, logg
    real,dimension(npatch) :: pcover
    integer,dimension(npatch):: do_clouds
    character(len=15),intent(inout) :: gasname(:)
    real, intent(inout) :: molmass(:)
    double precision, intent(inout) :: logVMR(:,:)
    character(len=15),INTENT(INOUT) :: cloudname(:,:)

    double precision,INTENT(INOUT) :: cloudrad(:,:,:)
    double precision,INTENT(INOUT) :: cloudsig(:,:,:)
    double precision,INTENT(INOUT) :: cloudprof(:,:,:)
    integer,intent(inout):: incloudnum(:,:)
    double precision,intent(inout) :: inwavenum(:)
    !real,dimension(nlinetemps) :: inlinetemps
    real,intent(inout) :: inlinetemps(:)
    real,intent(in) :: inpress(:)
    double precision,intent(inout):: linelist(:,:,:,:)
    double precision,intent(inout):: bff(:,:)
    real, dimension(nciatemps), intent(in):: ciatemp
    real, intent(inout) :: cia(:,:,:)
    double precision,allocatable, dimension(:,:),INTENT(OUT) :: out_spec
    double precision,allocatable, dimension(:,:),INTENT(INOUT) :: clphotspec, othphotspec
    double precision,allocatable, dimension(:,:,:),INTENT(INOUT) :: cf
    double precision,allocatable, dimension(:)::specflux
    real:: metal,grav,test,tau1
    ! counters
    integer :: ch4index,ipatch,icloud,ilayer,iwave,igas,nw1,nw2,use_disort
    integer :: do_bff
    real:: totcover, fboth, fratio, tstart,tfinish, opstart,opfinish,allelse
    real:: linstart, linfinish,distart, difinish,cloudstart,cloudfinish
    real:: bfstart,bffinish
    logical :: disorting,clphot,othphot,bfing,do_cf

    ! Are we using DISORT
    disorting = use_disort
    bfing  = do_bff
    ! HARD CODED to do line and CIA opacity calcs and everything
    ! apart from dust for first patch
    ! and copy to rest of patches before disort

    call init_all
    !wavenum[nwave],wavelen[nwave])


    do ipatch = 1, npatch
       call init_column(patch(ipatch)%atm)
    end do
    
    patch(1)%atm%temp = temp
   
    wavenum = inwavenum
    linetemps = inlinetemps
    press = inpress
    cloudnum = incloudnum

    call set_pressure_scale

    ! TEST LINE - set artificial temp profile
    !patch(1)%atm%temp = 300. + (((patch(1)%atm%press/1000.) / 0.1)**1.1)


    
    ch4index=0
    do igas = 1, ngas
       do ilayer = 1, nlayers
          
          patch(1)%atm(ilayer)%gas(igas)%name = trim(gasname(igas))
          patch(1)%atm(ilayer)%gas(igas)%VMR = 10.**(logVMR(igas,ilayer))
          patch(1)%atm(ilayer)%gas(igas)%molmass = molmass(igas)
       end do
       if ((trim(patch(1)%atm(1)%gas(igas)%name) .eq. "ch4") .or. (trim(patch(1)%atm(1)%gas(igas)%name) .eq. "12C1H4")) then
          ch4index = igas
       end if

    end do

    if (bfing) then
       patch(1)%atm%fe = 10.**bff(1,:)
       patch(1)%atm%fH = 10.**bff(2,:)
       patch(1)%atm%fHmin = 10.**bff(3,:)
    else
       patch(1)%atm%fe = 0.
       patch(1)%atm%fH = 0.
       patch(1)%atm%fHmin = 0.
    endif
   
    grav = 10**(logg) / 100.
    
    
    ! now let's check how long this takes
    
    call cpu_time(tstart)
    
    
    
    ! now H2 and He fractions and mu for each layer
    do ilayer = 1, nlayers

       allelse = sum(patch(1)%atm(ilayer)%gas%VMR) + patch(1)%atm(ilayer)%fe &
            + patch(1)%atm(ilayer)%fH + patch(1)%atm(ilayer)%fHmin
       !write(*,*) allelse
       fboth = 1.0 - allelse
       
       ! hardcoded H/He ratio
       ! from solar abundance of 91.2 by number H, 8.7 by number He
       fratio = 0.84
       patch(1)%atm(ilayer)%fH2 = fratio * fboth
       patch(1)%atm(ilayer)%fHe = (1.0  - fratio) * fboth
       
       patch(1)%atm(ilayer)%mu = (patch(1)%atm(ilayer)%fH2 * XH2) + &
            (patch(1)%atm(ilayer)%fHe * XHe) + &
            (patch(1)%atm(ilayer)%fH * XH) + &
            (patch(1)%atm(ilayer)%fHmin * XH) + &
            sum(patch(1)%atm(ilayer)%gas%VMR * patch(1)%atm(ilayer)%gas%molmass)
       !write(*,*) patch(1)%atm(ilayer)%mu
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
          patch(ipatch)%atm(ilayer)%opd_scat = 0d0
          patch(ipatch)%atm(ilayer)%gg = 0d0
          patch(ipatch)%atm(ilayer)%opd_CIA = 0d0
          patch(ipatch)%atm(ilayer)%opd_ext = 0d0
          patch(ipatch)%atm(ilayer)%opd_lines = 0d0
          patch(ipatch)%atm(ilayer)%opd_rayl = 0d0
          patch(ipatch)%atm(ilayer)%opd_hmbff = 0d0          
       end do
    end do
    
    wavelen = 1e4 / wavenum
    
     
    call cpu_time(opstart)
    ! now mix the gases in each layer to get optical depth from lines
    do ilayer = 1, nlayers
       call line_mixer(patch(1)%atm(ilayer),patch(1)%atm(ilayer)%opd_lines,ilayer,linelist)
    end do
    
    ! get the bff opacities
    if (bfing) then
       call get_hmbff
    end if
    
    ! now the Rayleigh scattering
    call get_ray(ch4index)

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
    call cpu_time(cloudstart)
    do ipatch = 1, npatch
       
       
       patch(ipatch)%cover = pcover(ipatch)
       
       patch(ipatch)%cloudy = do_clouds(ipatch)
       
       if (patch(ipatch)%cloudy .ne. 0) then
          do icloud = 1, nclouds
             patch(ipatch)%atm(1)%cloud(icloud)%name = cloudname(ipatch,icloud)
             
             
             ! in case of simple/generic/mixed cloud we won't be doing Mie coeffs
             ! we'll just use  rg, and rsig as w0 and gg
             ! for the cloud
             if (cloudnum(ipatch,icloud) .gt. 50) then
                !if (nclouds .ne. 1) then
                !   write(*,*) "Error: mixto cloud case should have nclouds = 1"
                !   stop
                !else
                ! FOR GREY SIMPLE CASE (MIXTO): put DTAU_cloud in cloudprofile
                ! and albedo in rg, and asymmetry = 0.0.
                ! Power law option has power law in rsig
                !
                do ilayer= 1, nlayers
                   patch(ipatch)%atm(ilayer)%opd_ext = &
                        patch(ipatch)%atm(ilayer)%opd_ext + &
                        ((cloudprof(ipatch,ilayer,icloud) * &
                        (wavelen**cloudsig(ipatch,ilayer,icloud))))
                   patch(ipatch)%atm(ilayer)%opd_scat = &
                        patch(ipatch)%atm(ilayer)%opd_scat + &
                        ((cloudprof(ipatch,ilayer,icloud) * &
                        (wavelen**cloudsig(ipatch,ilayer,icloud))* &
                         cloudrad(ipatch,ilayer,icloud)))
                   patch(ipatch)%atm(ilayer)%gg = 0.d0 
                end do ! layer loop
             else !if cloud < 50
                do ilayer = 1, nlayers
                   patch(ipatch)%atm(ilayer)%cloud(icloud)%name = &
                        patch(ipatch)%atm(1)%cloud(icloud)%name
                   patch(ipatch)%atm(ilayer)%cloud(icloud)%dtau1 = &
                        cloudprof(ipatch,ilayer,icloud)
                   patch(ipatch)%atm(ilayer)%cloud(icloud)%rg = cloudrad(ipatch,ilayer,icloud) 
                   patch(ipatch)%atm(ilayer)%cloud(icloud)%rsig = cloudsig(ipatch,ilayer,icloud)
                end do
             endif
          end do ! cloud loop
          if (any(cloudnum(ipatch,:) .lt. 50)) then
             !write(*,*) "calling cloud atlas"
             call cloudatlas(patch(ipatch)%atm,patch(ipatch)%cloudy)
          end if
          !do ilayer = 1,nlayers
          !   tau1 = tau1  + patch(1)%atm(ilayer)%cloud(1)%dtau1
          !end do
          !write(*,*)'cloud 1 total optical depth = ', tau1
          
       else
          do ilayer = 1, nlayers
             patch(ipatch)%atm(ilayer)%gg = 0.0
             patch(ipatch)%atm(ilayer)%opd_scat = 0.0
             patch(ipatch)%atm(ilayer)%opd_ext = 0.0
          end do
       end if
       
    end do ! patch loop
    call cpu_time(cloudfinish)
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

    
    call run_RT(specflux,clphotspec,othphotspec,cf,disorting,clphot,othphot,do_cf)

    allocate(out_spec(2,nwave))
    out_spec(1,:) = wavelen
    out_spec(2,:) = specflux * R2D2
    ! scale by r2d2
    
    call cpu_time(difinish)
    
    
    call cpu_time(tfinish)
    
    !write(*,*) "Time elapsed :", (tfinish - tstart), " seconds"
    
    !write(*,*) "Opacity interpolations took : ", (opfinish - opstart), " seconds"

    !write(*,*) "Cloud bits took: ", (cloudfinish - cloudstart), " seconds"
    !write(*,*) "RT took : ", (difinish - distart), " seconds"

    deallocate(wavelen,wavenum)
    do ipatch = 1, npatch
       do ilayer= 1, nlayers
          deallocate(patch(ipatch)%atm(ilayer)%opd_ext,&
               patch(ipatch)%atm(ilayer)%opd_scat,&
               patch(ipatch)%atm(ilayer)%opd_lines, &
               patch(ipatch)%atm(ilayer)%opd_cia, &
               patch(ipatch)%atm(ilayer)%opd_hmbff,&
               patch(ipatch)%atm(ilayer)%gg,&
               patch(ipatch)%atm(ilayer)%cloud,&
               patch(ipatch)%atm(ilayer)%gas)
       end do
       deallocate(patch(ipatch)%atm)
    end do
   
  end subroutine forward

  
end module main
