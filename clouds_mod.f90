module clouds

contains
  
  subroutine cloudatlas(column)

    use sizes
    use common_arrays
    use define_types
    use phys_const
   
    
    implicit none
    type(a_layer):: column(nlayers)
    integer :: icloud, imiewave, irad,ilayer,oldw1, oldw2, idum1, idum2, iwave
    character(len=50) :: miefile
    double precision, dimension(nmiewave,nrad,nclouds):: qscat,qext,cos_qscat
    double precision, dimension(nmiewave):: miewavelen,miewaven, wdiff
    double precision, dimension(nrad,nclouds) :: radius_in, radius, dr, rup
    double precision, dimension(nlayers,nmiewave,nclouds):: scat_cloud,ext_cloud,cqs_cloud
    double precision, dimension(nlayers,nmiewave):: opd_ext,opd_scat, cos_qs
    double precision :: norm, rr, r2, rg, rsig, pw, pir2ndz, arg1, arg2
    double precision :: f1, f2, intfact,lintfact


    ! this will take the clouds and in turn calculate their opacities the layers
    ! based on name, density, mean radius and width of log normal distribution


    ! first set up the grids and get the Mie coefficients, cloud by cloud

    call init_column(column)
    do icloud = 1, nclouds
       
       ! This bit is lifted from setup_clouds3.1.f in EGP
       ! it sets up the radius grid.
       
       ! constants/ parameters are in sizes
       
       pw = 1. / 3.
       f1 = ( 2*vrat / ( 1 + vrat) )**pw
       f2 = ( 2 / ( 1 + vrat ) )**pw * vrat**(pw-1)
       
       do irad = 1, nrad
          radius(irad,icloud) = rmin * vrat**(float(irad-1)/3.)
          rup(irad,icloud) = f1*radius(irad,icloud)
          dr(irad,icloud) = f2*radius(irad,icloud)
       enddo
       
       ! first get the mie coefficients for the condensate
       
       write(miefile,"(A,A,A)")"../Clouds/",trim(column(1)%cloud(icloud)%name),".mieff"
       
       
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
    end do ! cloud loop
    
    miewaven = 1.0 / miewavelen
    
    
    ! now for each layer, run through the particle sizes summing up the
    ! cross sections to get optical depth for this cloud from extinction
    ! and scattering
    ! this is all hacked from calc_optics. Credit to Ackerman & Marley
    scat_cloud = 0.d0
    ext_cloud = 0.d0
    cqs_cloud = 0.d0
    opd_ext = 0.d0
    opd_scat = 0.0
    cos_qscat = 1.d-99
    
    do ilayer =1, nlayers
       do icloud = 1, nclouds
          if (column(ilayer)%cloud(icloud)%density .gt. 0.) then
             
             
             rsig = column(ilayer)%cloud(icloud)%rsig
             rg  = column(ilayer)%cloud(icloud)%rg
             
             r2 = rg**2 * exp( 2*log(rsig)**2 )
             
             ! check the logic for this bit!!!
             ! Optical depth for conservative geometric scatterers 
             !opd_layer(ilayer,icloud) = 2.*pi*r2 * &
             !     column(ilayer)%cloud(icloud)%density * &
             !     column(ilayer)%dz
             
             !  Calculate normalization factor (forces lognormal sum = 1.0)
             
             norm = 0.
             
             do irad = 1,nrad
                rr = radius(irad,icloud)
                arg1 = dr(irad,icloud) / ( sqrt(2.*PI)*rr*log(rsig) )
                arg2 = -log( rr/ rg )**2 / ( 2*log(rsig)**2 )
                norm = norm + arg1*exp( arg2 )
             end do
             
             ! my lengths and densities are in metres, but radii are in cm
             norm = (column(ilayer)%cloud(icloud)%density * column(ilayer)%dz  * 10.**(-4)) / norm
             
             
             ! now loop over radius and fill up wavelength dependent opacity for
             ! each cloud
             do imiewave = 1, nmiewave
                do irad = 1, nrad
                   
                   rr = radius(irad,icloud)
                   arg1 = dr(irad,icloud) / ( sqrt(2.*PI)*log(rsig) )
                   arg2 = -log( rr/rg)**2 / ( 2*log(rsig)**2 )
                   pir2ndz = norm*PI*rr*arg1*exp( arg2 )
                                      
                                      
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
                     (cqs_cloud(ilayer,imiewave,icloud) * &
                     scat_cloud(ilayer,imiewave,icloud))
      
             end do ! miewave loop
          end if
       end do  ! cloud loop

       cos_qs(ilayer,:) = cos_qs(ilayer,:) / opd_scat(ilayer,:)

 
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
          
          intfact = (log10(wavenum(iwave)) - log10(miewaven(oldw1))) / &
               (log10(miewaven(oldw2)) - log10(miewaven(oldw1)))

          ! gg should be interpolated in linear space - it is always small
          lintfact =  (wavenum(iwave) - miewaven(oldw1)) / &
               (miewaven(oldw2) - miewaven(oldw1))
          
          column(ilayer)%opd_ext(iwave) = 10.**&
               (((log10(opd_ext(ilayer,oldw2)) - log10(opd_ext(ilayer,oldw1))) *intfact) &
               + log10(opd_ext(ilayer,oldw1)))
          
          column(ilayer)%opd_scat(iwave) = 10.**&
               (((log10(opd_scat(ilayer,oldw2)) - log10(opd_scat(ilayer,oldw1))) *intfact) &
               + log10(opd_scat(ilayer,oldw1)))
          
          column(ilayer)%gg(iwave) = &
               ((cos_qs(ilayer,oldw2) - cos_qs(ilayer,oldw1)) *lintfact) &
               + cos_qs(ilayer,oldw1) 
          
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
! write(*,*) "clouds line 187 opd_scat layer ",column(1:20)%opd_scat(300)
! write(*,*) "clouds line 188 opd_ext layer ",column(1:20)%opd_ext(300)
! write(*,*) "clouds line 189 gg layer ",column(1:20)%gg(300)

 
 
end subroutine cloudatlas



  
end module clouds
