subroutine properties(press,wavenum,nclouds,do_clouds,cloudnum,cloudprof,cloudrad,cloudsig,tau1_cloud,mass_cloud,num_cloud)




  !f2py intent(in) :: nclouds
  !f2py intent(inout) :: wavenum,press, do_clouds, cloudnum
  !f2py intent(inout) :: cloudprof,cloudrad,cloudsig
  !f2ly intent(out) :: mass_cloud, num_cloud, tau1_cloud
  
  implicit none
  integer,parameter :: nmiewave = 196,nrad = 60
  integer :: icloud, imiewave, irad,ilayer,oldw1, oldw2, idum1, idum2,iwave
  integer :: sizdist,loc1,loc1a,loc1b,loc1m,nwave,npatch,ipatch,nlayers,maxcloud
  integer :: loc(1)
  character(len=50) :: miefile
  character(len=15), allocatable, dimension(:,:) :: cloudname
  character(len=15),allocatable,dimension(:) :: cloudlist
  double precision, allocatable,dimension(:) :: dp
  double precision ,allocatable, dimension(:,:,:) :: qscat,qext,cos_qscat
  double precision ,allocatable, dimension(:,:,:,:) :: ext_cloud,tau_cloud
  double precision, dimension(196):: miewavelen,miewaven, wdiff
  double precision, allocatable, dimension(:,:) :: radius_in, radius
  double precision, allocatable, dimension(:,:) :: dr, rup
  logical, allocatable, dimension(:,:) :: cldone 
  double precision :: norm, rr, r2, rg, rsig, pw, pir2ndz, arg1, arg2,bot
  double precision :: f1, f2, intfact,lintfact,vrat,rmin, PI
  double precision :: a, b, ndz, drr, arg3, argnum, argext, argmass
  double precision :: logcon, qpir2, frac,rho,numdens,p1,p2,tau1,tau2
  integer,intent(in) :: nclouds
  double precision,intent(inout) :: press(:), wavenum(:)
  integer,intent(inout) :: do_clouds(:)
  integer,intent(inout) :: cloudnum(:,:)
  double precision,intent(inout) :: cloudrad(:,:,:)
  double precision,intent(inout) :: cloudsig(:,:,:)
  double precision,intent(inout) :: cloudprof(:,:,:)
  ! declare cloud mass and number arrays as max 2 patch,100 layers,6 clouds
  double precision,dimension(2,100,6),intent(OUT):: mass_cloud,num_cloud
  ! declare tau=1 array as 2 patch, 60000 wave points 6 clouds 
  double precision,dimension(2,60000,6),intent(OUT):: tau1_cloud
  
  PI = 3.14159274d0
  

  
  ! this will take the clouds and in turn calculate their opacities the layers
  ! based on name, density, mean radius and width of log normal distribution
  ! or using Hansen distribution: effective radius, radius spread and density

  npatch = size(do_clouds)
  nlayers = size(press)
  nwave = size(wavenum)
  
 
  
  allocate(qscat(nmiewave,nrad,nclouds),qext(nmiewave,nrad,nclouds))
  allocate(cos_qscat(nmiewave,nrad,nclouds))
  allocate(radius_in(nrad,nclouds), radius(nrad,nclouds))
  allocate(dr(nrad,nclouds), rup(nrad,nclouds))
  allocate(ext_cloud(npatch,nlayers,nmiewave,nclouds),tau_cloud(npatch,nlayers,nwave,nclouds))
  allocate(cloudname(npatch,nclouds))
  allocate(cldone(nclouds,nwave))
  allocate(dp(nlayers))

  open(10,file="data/cloudlist.dat", status='old')
  read(10,*) maxcloud
  allocate(cloudlist(maxcloud))
  do icloud = 1, maxcloud
     read(10,"(I3,A15)") idum2, cloudlist(icloud)
  end do
  close(10)

   
  
  do ipatch = 1, npatch
     do icloud = 1, nclouds
        ! check if we're doing a specific cloud or a generic/mixed cloud
        if (cloudnum(ipatch,icloud) .gt. 50) then
           cloudname(ipatch,icloud) = "mixto"
        else
           cloudname(ipatch,icloud) = trim(adjustl(cloudlist(cloudnum(ipatch,icloud))))
        endif
     end do
  end do

  ! first set up the grids and get the Mie coefficients, cloud by cloud
  
  do ipatch = 1, npatch
     cldone = .false.
     
     ! get the distribution from do_clouds
     ! 2 = log normal, 1 = hansen
     sizdist = do_clouds(ipatch)
     do icloud = 1, nclouds
     
        ! This bit is lifted from setup_clouds3.1.f in EGP
        ! it sets up the radius grid.
     
        ! constants/ parameters are in sizes
        
        ! set rmin and vrat here though:
     
        if (trim(cloudname(ipatch,icloud)) .eq. 'soot' ) then
           vrat = 1.3
           rmin = 1e-8
        else
           vrat = 2.2
           rmin = 1e-7
        endif
     
     
        pw = 1. / 3.
        f1 = ( 2*vrat / ( 1 + vrat) )**pw
        f2 = ( 2 / ( 1 + vrat ) )**pw * vrat**(pw-1)
     
        do irad = 1, nrad
           radius(irad,icloud) = rmin * vrat**(float(irad-1)/3.)
           rup(irad,icloud) = f1*radius(irad,icloud)
           dr(irad,icloud) = f2*radius(irad,icloud)
        enddo
     
        ! first get the mie or DHS  coefficients for the condensate
     
        if (sizdist .gt. 0) then
           write(miefile,"(A,A,A)")"../Clouds/",trim(cloudname(ipatch,icloud)),".mieff"
        else if (sizdist .lt. 0) then
           write(miefile,"(A,A,A)")"../Clouds/",trim(cloudname(ipatch,icloud)),".dhs"
        end if
     
     
        open(unit=10, file= miefile, status='old')
     
     
        read(10,*) idum1, idum2
     
        if (idum1 .ne. nmiewave .or. idum2 .ne. nrad) then
           write(*,*) "Problem with mie coefficients file contents wrong waves or radii",trim(miefile)
           stop
        end if
     
        do irad = 1, nrad
           read(10,*) radius_in(irad,icloud)
           if (abs(radius(irad,icloud) - radius_in(irad,icloud)) .gt. &
                0.01*radius(irad,icloud)) then
              write(*,*) "Radius grid mismatch in mie file: ",trim(miefile)
              stop
           end if
           do imiewave = 1, nmiewave
              read(10,*) miewavelen(imiewave), qscat(imiewave,irad,icloud), &
                   qext(imiewave,irad,icloud), cos_qscat(imiewave,irad,icloud)
           end do
        end do
        close(10)   
     end do  ! cloud loop

  
     miewaven = 1.0 / miewavelen
  
     ! let's get the location for wave = 1um in miewavelen for later
     loc = minloc(abs(miewavelen - 1e-4))
     loc1 = loc(1)
     !write(*,*) miewavelen(loc1)
  
  
  
     ! now for each layer, run through the particle sizes summing up the
     ! cross sections to get optical depth for this cloud from extinction
     ! and scattering
     ! this is all hacked from calc_optics. Credit to Ackerman & Marley
     
     ! ensure that we're starting from zero values 
     ext_cloud = 0.d0
     mass_cloud = 0.d0
     num_cloud = 0.d0


     do ilayer = 1, nlayers
        ! layer dP
        if (ilayer .eq. nlayers) then
          p1 =  exp((0.5)*(log(press(ilayer-1) * press(ilayer))))

          p2 = press(ilayer)**2 / p1

          dp(ilayer) = p2 - p1
       else
          p1 = exp(((1.5)*log(press(ilayer)) - &
               ((0.5)*log(press(ilayer+1)))))
          p2 = exp((0.5)*(log(press(ilayer) * press(ilayer+1))))
          
          dp(ilayer) = p2 - p1
           
       end if

       do icloud = 1, nclouds
          ! first we need the cloud density in this layer
          ! we get this from the layer optical depth of the cloud at 1um
          ! which is what we're given
          
          if (cloudprof(ipatch,ilayer,icloud) .gt. 1.d-6) then
             if (abs(sizdist) .eq. 2) then
                ! we take geometric mean parameter from python code
                ! as a value between 0 and 1. This is then translated here to
                ! hold a value between 1 and 5
                rsig = 1. + (cloudsig(ipatch,ilayer,icloud) * 4)
                ! radii supplied in um, convert to cm
                rg  = cloudrad(ipatch,ilayer,icloud) * 1e-4
                
                r2 = rg**2 * exp( 2*log(rsig)**2 )
                
                !  Calculate normalization factor , i.e N_0
                ! This is based on setting tau_cl = 1 at 1 micron
                ! so we sum up the cross-section contrbutions at 1um from
                ! all radii particles across the distribution
                ! get Ndz from 1/ this sum
                norm = 0.
                
                do irad = 1,nrad
                   rr = radius(irad,icloud)
                   arg1 = dr(irad,icloud) / ( sqrt(2.*PI)*rr*log(rsig) )
                   arg2 = -log( rr/ rg )**2 / ( 2*log(rsig)**2 )
                   qpir2 = PI * rr**2 * qext(loc1,irad,icloud)
                   norm = norm + (qpir2 * arg1 * exp(arg2))
                end do
                
                ! so Ndz (i.e total number density * height of layer) 
                ndz  =  1. / norm
                
                
                ! now loop over radius and fill up wavelength dependent opacity for
                ! each cloud
                do imiewave = 1, nmiewave
                   do irad = 1, nrad
                      
                      rr = radius(irad,icloud)
                      arg1 = dr(irad,icloud) / ( sqrt(2.*PI)*rr*log(rsig) )
                      arg2 = -log( rr/rg)**2 / ( 2*log(rsig)**2 )
                      pir2ndz = ndz * PI * rr**2 * arg1* exp( arg2 )
                      
                      
                      ext_cloud(ipatch,ilayer,imiewave,icloud) = &
                           ext_cloud(ipatch,ilayer,imiewave,icloud) + &
                           qext(imiewave,irad,icloud)*pir2ndz
                   enddo ! radius loop                
                end do ! miewave loop
                
             elseif (abs(sizdist) .eq. 1) then
                
                ! Hansen distribution
                
                ! radii supplied in um, convert to cm
                a = cloudrad(ipatch,ilayer,icloud) * 1d-4
                ! b is not a length, it is dimensionless
                b  = cloudsig(ipatch,ilayer,icloud)
              
                ! first need to get ndz from the optical depth dtau at 1um
                 
                bot = 0.d0
                do irad = 1, nrad
                   rr = radius(irad,icloud)
                   drr = dr(irad,icloud)
                   !write(*,*) (-rr/(a*b)), log(drr)
                   arg1 = (-rr/(a*b)) + log(drr)
                   !write(*,*) arg1
                   arg2 = ((1.- 3.*b)/b) * log(rr)
                   argext = log(qext(loc1,irad,icloud) * PI * rr**2.)
                   bot = bot + exp(arg1 + arg2 + argext)
                   
                end do ! radius loop
                
                logcon = log(cloudprof(ipatch,ilayer,icloud) / bot)
                
                arg3 = ((((2.*b) - 1.)/b) * log(a*b))
                arg2 = log_gamma((1.-(2.*b))/b)
                
                
                ndz = exp(logcon +arg2 - arg3)
                
                arg1 = ((((2.*b) - 1.)/b) * log(a*b)) + log(ndz)
                
                logcon =  (arg1 - arg2) 

                do imiewave = 1, nmiewave
                   do irad = 1, nrad
                      rr = radius(irad,icloud)
                      drr = dr(irad,icloud)
                      !write(*,*) (-rr/(a*b)), log(drr)
                      arg1 = (-rr/(a*b)) + log(drr)
                      !write(*,*) arg1
                      arg2 = ((1. - 3.*b)/b) * log(rr)
                      argext = log(qext(imiewave,irad,icloud) * PI * rr**2)
                      !write(*,*) logcon, arg1, arg2, arg3, arg4 
                      
                      ext_cloud(ipatch,ilayer,imiewave,icloud) = &
                           ext_cloud(ipatch,ilayer,imiewave,icloud) + &
                           exp(logcon + arg1 + arg2 + argext)
                      
                      
                      ! now let's get masses
                      ! just for MgSiO3 and SiO2 for now....
                      ! just do for miewave == 1. 
                      
                      if (imiewave .eq. 1) then
                         if (cloudnum(ipatch,icloud) .eq. 5) then
                            rho = 3.66 ! g /cm3
                            ! https://onlinelibrary.wiley.com/doi/full/10.1111/j.1945-5100.2010.01129.x
                            
                            numdens = 2.196e22 ! "molecules per cm3"
                            
                         elseif (cloudnum(ipatch,icloud) .eq. 12) then
                            rho = 2.65 ! g/cm3 from google
                            
                            numdens = 2.6560e22 ! molecules per cm3
                         end if
                         
                         argmass = log(rho * (4/3) * PI * rr**3)
                         argnum = log(numdens * (4/3) * PI * rr**3) 
                         mass_cloud(ipatch,ilayer,icloud) = &
                              mass_cloud(ipatch,ilayer,icloud) + &
                              exp(logcon + arg1 + arg2 + argmass)
                         
                         num_cloud(ipatch,ilayer,icloud) = &
                              num_cloud(ipatch,ilayer,icloud) + &
                              exp(logcon + arg1 + arg2 + argnum)
                      end if
                   end do  ! radius loop
                   
                end do ! miewave loop
             end if
          end if
       end do   ! cloud loop
     
       ! rebin to working resolution (nwave) grid 
       do iwave= 1 , nwave




          wdiff = abs(miewaven - wavenum(iwave))
          
          oldw1 = minloc(wdiff,1)
          
          if (miewaven(oldw1) .lt. wavenum(iwave)) then
             oldw2 = oldw1 + 1
          else
             oldw2 = oldw1
             oldw1 = oldw2 - 1
          end if

          !intfact = (log10(wavenum(iwave)) - log10(miewaven(oldw1))) / &
          !     (log10(miewaven(oldw2)) - log10(miewaven(oldw1)))
          
          ! gg should be interpolated in linear space - it is always small
          lintfact =  (wavenum(iwave) - miewaven(oldw1)) / &
               (miewaven(oldw2) - miewaven(oldw1))
          
          !loop through clouds here
          do icloud = 1, nclouds
             tau_cloud(ipatch,ilayer,iwave,icloud) = &
                  ((ext_cloud(ipatch,ilayer,oldw2,icloud) - ext_cloud(ipatch,ilayer,oldw1,icloud)) * lintfact) &
                  + ext_cloud(ipatch,ilayer,oldw1,icloud)
             
             if (tau_cloud(ipatch,ilayer,iwave,icloud) .lt. 1d-50) then
                tau_cloud(ipatch,ilayer,iwave,icloud) = 0.
             end if

             ! now we;re getting the tau = 1 pressure as a function of wavenum for each cloud
             if (cldone(icloud,iwave) .eqv. .false.) then
                if (sum(tau_cloud(ipatch,1:ilayer,iwave,icloud)) .gt. 1.0) then
                   cldone(icloud,iwave) = .true.
                   tau2 = sum(tau_cloud(ipatch,1:ilayer,iwave,icloud))
                   tau1 = tau2 - tau_cloud(ipatch,ilayer,iwave,icloud)
                   if (ilayer .eq. nlayers) then
                      p1 = exp((0.5)*(log(press(ilayer-1) * press(ilayer))))
                   else
                      p1 = exp(((1.5)*log(press(ilayer))) - &
                        ((0.5)*log(press(ilayer+1))))
                   end if
                   tau1_cloud(ipatch,iwave,icloud) = p1 +((1.0 - tau1) * &
                     dp(ilayer) / tau_cloud(ipatch,ilayer,iwave,icloud))
                end if
             end if
          end do ! cloud loop
       end do ! wave loop
    end do  ! layer loop
 end do ! patch loop


  
  
 deallocate(cloudlist)
 deallocate(qscat,qext)
 deallocate(cos_qscat)
 deallocate(radius_in,radius)
 deallocate(dr,rup)
 deallocate(ext_cloud,tau_cloud)
 deallocate(cloudname)
 deallocate(cldone)
 deallocate(dp)

 
  
  
end subroutine properties
