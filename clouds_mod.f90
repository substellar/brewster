module clouds

contains
  
  subroutine cloudatlas(column,sizdist)

    use sizes
    use common_arrays
    use define_types
    use phys_const
   
    
    implicit none
    type(a_layer), intent(inout):: column(nlayers)
    integer :: icloud, imiewave, irad,ilayer,oldw1, oldw2, idum1, idum2,iwave
    integer :: sizdist,loc1,loc1a,loc1b,loc1m
    integer :: loc(1)
    character(len=50) :: miefile
    double precision ,allocatable, dimension(:,:,:) :: qscat,qext,cos_qscat
    double precision, dimension(nmiewave):: miewavelen,miewaven, wdiff,ext,scat,cqs
    double precision, allocatable, dimension(:,:) :: radius_in, radius
    double precision, allocatable, dimension(:,:) :: dr, rup
    double precision, allocatable, dimension(:,:,:) :: scat_cloud,ext_cloud
    double precision, allocatable, dimension(:,:,:) :: cqs_cloud
    double precision, allocatable, dimension(:,:):: opd_ext,opd_scat, cos_qs
    double precision :: norm, rr, r2, rg, rsig, pw, pir2ndz, arg1, arg2,bot
    double precision :: f1, f2, intfact,lintfact,vrat,rmin,sum_ndz_2, sum_ndz_1
    double precision :: a, b, ndz, drr, arg3, argscat, argext, argcosqs
    double precision :: logcon, qpir2, frac
    real, allocatable,dimension(:) :: cld1arr 
    ! this will take the clouds and in turn calculate their opacities the layers
    ! based on name, density, mean radius and width of log normal distribution
    ! or using Hansen distribution: effective radius, radius spread and density

    ! get the distribution from do_clouds
    ! 2 = log normal, 1 = hansen
    allocate(qscat(nmiewave,nrad,nclouds),qext(nmiewave,nrad,nclouds))
    allocate(cos_qscat(nmiewave,nrad,nclouds))
    allocate(radius_in(nrad,nclouds), radius(nrad,nclouds))
    allocate(dr(nrad,nclouds), rup(nrad,nclouds))
    allocate(scat_cloud(nlayers,nmiewave,nclouds))
    allocate(ext_cloud(nlayers,nmiewave,nclouds))
    allocate(cqs_cloud(nlayers,nmiewave,nclouds))
    allocate(opd_ext(nlayers,nmiewave),opd_scat(nlayers,nmiewave))
    allocate(cos_qs(nlayers,nmiewave))
    allocate(cld1arr(nlayers))

    ! first set up the grids and get the Mie coefficients, cloud by cloud

    !call init_column(column)
    do icloud = 1, nclouds
       
       ! This bit is lifted from setup_clouds3.1.f in EGP
       ! it sets up the radius grid.
       
       ! constants/ parameters are in sizes
       
       ! set rmin and vrat here though:
       
       if (trim(column(1)%cloud(icloud)%name) .eq. 'soot' ) then
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
          write(miefile,"(A,A,A)")"../Clouds/",trim(column(1)%cloud(icloud)%name),".mieff"
       else if (sizdist .lt. 0) then
          write(miefile,"(A,A,A)")"../Clouds/",trim(column(1)%cloud(icloud)%name),".dhs"
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
    scat_cloud = 0.d0
    ext_cloud = 0.d0
    cqs_cloud = 0.d0
    opd_ext = 0.d0
    opd_scat = 0.d0
    cos_qs = 0.d0
    sum_ndz_1 = 0.d0
    sum_ndz_2 = 0.d0
    ! This sets up options to print ndz etc
    do ilayer = 1, nlayers
       cld1arr(ilayer) = column(ilayer)%cloud(1)%dtau1
    end do
    loc = maxloc(cld1arr)
    idum1 = loc(1)
 

    do ilayer = 1, nlayers
       do icloud = 1, nclouds
          ! first we need the cloud density in this layer
          ! we get this from the layer optical depth of the cloud at 1um
          ! which is what we're given
          
          if (column(ilayer)%cloud(icloud)%dtau1 .gt. 1.d-6) then

             if (abs(sizdist) .eq. 2) then
                ! we take geometric mean parameter from python code
                ! as a value between 0 and 1. This is then translated here to
                ! hold a value between 1 and 5
                rsig = 1. + (column(ilayer)%cloud(icloud)%rsig * 4)
                ! radii supplied in um, convert to cm
                rg  = column(ilayer)%cloud(icloud)%rg * 1e-4
                
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
                   arg2 = -(log( rr/ rg ))**2 / ( 2*(log(rsig))**2 )
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
                      arg2 = -(log( rr/rg))**2 / ( 2*(log(rsig))**2 )
                      pir2ndz = ndz * PI * rr**2 * arg1* exp( arg2 )
                      
                      
                      scat_cloud(ilayer,imiewave,icloud) =  &
                           scat_cloud(ilayer,imiewave,icloud) + & 
                           qscat(imiewave,irad,icloud)*pir2ndz
                      ext_cloud(ilayer,imiewave,icloud) = &
                           ext_cloud(ilayer,imiewave,icloud) + &
                           qext(imiewave,irad,icloud)*pir2ndz
                      cqs_cloud(ilayer,imiewave,icloud) = &
                           cqs_cloud(ilayer,imiewave,icloud) + &
                           cos_qscat(imiewave,irad,icloud)*pir2ndz
                   enddo ! radius loop
                   
                   ! sum over clouds
                   opd_scat(ilayer,imiewave) = opd_scat(ilayer,imiewave) + &
                        scat_cloud(ilayer,imiewave,icloud)
                   opd_ext(ilayer,imiewave) = opd_ext(ilayer,imiewave) + &
                        ext_cloud(ilayer,imiewave,icloud)
                   cos_qs(ilayer,imiewave) = cos_qs(ilayer,imiewave) + &
                        cqs_cloud(ilayer,imiewave,icloud)
!                        scat_cloud(ilayer,imiewave,icloud))
      
                end do ! miewave loop
             elseif (abs(sizdist) .eq. 1) then

                ! Hansen distribution
                
                ! radii supplied in um, convert to cm
                a = column(ilayer)%cloud(icloud)%rg * 1d-4
                ! b is not a length, it is dimensionless
                b  = column(ilayer)%cloud(icloud)%rsig
                

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
                
                logcon = log(column(ilayer)%cloud(icloud)%dtau1 / bot)

                arg3 = ((((2.*b) - 1.)/b) * log(a*b))
                arg2 = log_gamma((1.-(2.*b))/b)

                
                ndz = exp(logcon +arg2 - arg3)


                !if (icloud .eq. 1) then
                !   sum_ndz_1 = sum_ndz_1 +ndz
                !   frac = (ndz *1e4 / (column(ilayer)%dz)) * ((4/3)*PI * (0.008e-4**3)) * 0.032 * AVOGADRO / (column(ilayer)%ndens)
                !   write(*,*) 'MgSiO3 mixing fraction = ', frac
                ! write(*,*) 'sum ndz MgSiO3 = ', sum_ndz_1
                !end if
                !if (icloud .eq. 2) then
                !   sum_ndz_2 = sum_ndz_2 +ndz
                !   frac = (ndz * 1e4/ (column(ilayer)%dz)) * ((4/3)*PI * (0.34e-4**3)) * 0.044 * AVOGADRO / (column(ilayer)%ndens)
                !   write(*,*) 'SiO2 mixing fraction = ', frac
                !   write(*,*) 'sum ndz SiO2 = ', sum_ndz_2
                !end if

                !TEST writes for ndz etc
                !if (ilayer .eq. idum1) then
                !   write(*,*) 'ndz = ', ndz
                !   write(*,*)  'dz = ', column(ilayer)%dz
                !   write(*,*) 'n = ',ndz / column(ilayer)%dz,' m^-3'
                !   write(*,*) 'at pressure ', column(ilayer)%press,' bars'
                !end if
                
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
                      argscat = log(qscat(imiewave,irad,icloud) * PI * rr**2)
                      argext = log(qext(imiewave,irad,icloud) * PI * rr**2)
                      argcosqs =  cos_qscat(imiewave,irad,icloud) * PI * rr**2
                      !write(*,*) logcon, arg1, arg2, arg3, arg4 
                      scat_cloud(ilayer,imiewave,icloud) =  &
                           scat_cloud(ilayer,imiewave,icloud) + &
                           exp(logcon + arg1 + arg2 + argscat)
                      
                      ext_cloud(ilayer,imiewave,icloud) = &
                           ext_cloud(ilayer,imiewave,icloud) + &
                           exp(logcon + arg1 + arg2 + argext)
                      
                      cqs_cloud(ilayer,imiewave,icloud) = &
                           cqs_cloud(ilayer,imiewave,icloud) + &
                           (exp(logcon + arg1 + arg2) * argcosqs)

                   end do  ! radius loop                                     
                   ! sum over clouds
                   opd_scat(ilayer,imiewave) = opd_scat(ilayer,imiewave) + &
                        scat_cloud(ilayer,imiewave,icloud)
                   opd_ext(ilayer,imiewave) = opd_ext(ilayer,imiewave) + &
                        ext_cloud(ilayer,imiewave,icloud)
                   cos_qs(ilayer,imiewave) = cos_qs(ilayer,imiewave) + &
                        cqs_cloud(ilayer,imiewave,icloud)
                end do ! miewave loop
             end if
             !write (*,*) scat_cloud(ilayer,loc1,icloud)
             !write (*,*) ext_cloud(ilayer,loc1,icloud)
             !write (*,*) cqs_cloud(ilayer,loc1,icloud)
             
          end if
       end do   ! cloud loop
       
       ! rebin to working resolution (nwave) grid and write to
       
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
          
          column(ilayer)%opd_ext(iwave) = &
               ((opd_ext(ilayer,oldw2) - opd_ext(ilayer,oldw1)) * lintfact) &
               + opd_ext(ilayer,oldw1)
          
          column(ilayer)%opd_scat(iwave) = &
               ((opd_scat(ilayer,oldw2) - opd_scat(ilayer,oldw1)) * lintfact) &
               + opd_scat(ilayer,oldw1)
          
          column(ilayer)%gg(iwave) = &
                  (((cos_qs(ilayer,oldw2) - cos_qs(ilayer,oldw1)) * lintfact) &
                  + cos_qs(ilayer,oldw1) ) / column(ilayer)%opd_scat(iwave)
                       
          if (column(ilayer)%opd_scat(iwave) .lt. 1d-50) then
             column(ilayer)%opd_scat(iwave) = 0.
             column(ilayer)%gg(iwave) = 0.
          end if
          if (column(ilayer)%opd_ext(iwave) .lt. 1d-50) then
             column(ilayer)%opd_ext(iwave) = 0.
          end if
          
          
       end do ! wave loop

    end do  ! layer loop

    ! TK test line see what we've got
    
    !loc = minloc(abs((1e4/wavenum) - 1.2 ))
    !loc1a = loc(1)
    
    !write(*,*) "wavenum, wavelength for check = ", wavenum(loc1a), wavelen(loc1a)
    !write(*,*) "clouds line 187 opd_scat layer ",column(idum1)%opd_scat(loc1a)
    !write(*,*) "clouds line 188 opd_ext layer",column(idum1)%opd_ext(loc1a)
    !write(*,*) "clouds line 189 gg layer ",column(idum1)%gg(loc1a)


    ! test line write cloud optical depth out
    !open(unit=10,file='cloud3_diagn.txt',form='formatted',status='new')
    !scat = 0.
    !ext = 0.
    !cqs = 0.
    !loc = minloc(abs((1e4*miewavelen) - 1.0 ))
    !loc1m = loc(1)
    !do ilayer = 1, nlayers
    !   if (ext(loc1m) .lt. 1.0) then
    !      do imiewave= 1, nmiewave
    !         scat(imiewave) = scat(imiewave) + scat_cloud(ilayer,imiewave,3)
    !         ext(imiewave) = ext(imiewave) + ext_cloud(ilayer,imiewave,3)
    !         cqs(imiewave) = cqs(imiewave) + cqs_cloud(ilayer,imiewave,3)
    !      end do
    !   end if
    !enddo
    !do imiewave = 1, nmiewave
    !   write(10,*) miewavelen(imiewave),scat(imiewave),ext(imiewave),cqs(imiewave)
    !end do
    !close(10)


    deallocate(qscat,qext,cos_qscat,radius_in,radius,dr,rup)
    deallocate(scat_cloud,ext_cloud,cqs_cloud)
    deallocate(opd_ext,opd_scat,cos_qs)
 
end subroutine cloudatlas



  
end module clouds
