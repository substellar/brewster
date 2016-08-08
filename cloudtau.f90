subroutine cloudtau(theta,cloudname,sizdist,miewavelen,tau)

  use sizes

  implicit none

  !f2py integer, parameter :: nmiewave

  !f2py intent(in) theta
  !f2py intent(in) cloudname
  !f2py intent(in) sizdist
  !f2py intent(out) tau
  !f2py intent(out) miewavelen

  integer :: imiewave, irad, idum1, idum2,sizdist
  character(len=50) :: miefile,cloudname
  double precision, dimension(nmiewave,nrad):: qscat,qext,cos_qscat
  double precision, dimension(nmiewave):: miewavelen
  double precision, dimension(nrad) :: radius_in, radius, dr, rup
  double precision, dimension(nmiewave):: ext_cloud,tau
  double precision,dimension(3)::theta
  double precision :: norm, rr, r2, rg, rsig, pw, pir2ndz, arg1, arg2,arg3
  double precision :: f1, f2,con,vrat,rmin,a,b,ndz,logcon,arg4,drr
  double precision, parameter :: pi = 3.14159274d0


  ! first set up the grids and get the Mie coefficients, cloud by cloud

  ! This bit is lifted from setup_clouds3.1.f in EGP
  ! it sets up the radius grid.
  

  
  ! constants/ parameters are in sizes

  if( trim(cloudname) .eq. 'CH4' )then
     vrat = 2.2
     rmin = 1e-5
  elseif( trim(cloudname) .eq. 'NH3' )then
     vrat = 2.2
     rmin = 1e-5
  elseif( trim(cloudname) .eq. 'H2O' )then
     vrat = 2.2
     rmin = 1e-5
  elseif( trim(cloudname) .eq. 'Fe' )then
     vrat = 2.2
     rmin = 1e-7
  elseif( trim(cloudname) .eq. 'KCl' )then
     vrat = 2.2
     rmin = 1e-5
  elseif( trim(cloudname) .eq. 'NaCl' )then
     vrat = 2.2
     rmin = 1e-5
  elseif( trim(cloudname) .eq. 'Cr' )then
     vrat = 2.2
     rmin = 1e-5
  elseif( trim(cloudname) .eq. 'MgSiO3' )then
     vrat = 2.2
     rmin = 1e-7
  elseif( trim(cloudname) .eq. 'Mg2SiO4' )then
     vrat = 2.2
     rmin = 1e-7
  elseif( trim(cloudname) .eq. 'ZnS' )then
     vrat = 2.2
     rmin = 1e-5
  elseif( trim(cloudname) .eq. 'MnS' )then
     vrat = 2.2
     rmin = 1e-7
  elseif( trim(cloudname) .eq. 'Na2S' )then
     vrat = 2.2
     rmin = 1e-7
  elseif( trim(cloudname) .eq. 'NH4H2PO4' )then
     vrat = 2.2
     rmin = 1e-5
  elseif( trim(cloudname) .eq. 'tholins' )then
     vrat = 2.2
     rmin = 1e-7
  elseif( trim(cloudname) .eq. 'soot' )then
     vrat = 1.3
     rmin = 1e-7
  elseif( trim(cloudname) .eq. 'testgrid3' )then
     vrat = 2.2
     rmin = 1e-5
  else
     write(*,*) 'init_optics(): bad igas = ', trim(cloudname)
     stop 
  endif

  
  pw = 1. / 3.
  f1 = ( 2*vrat / ( 1 + vrat) )**pw
  f2 = ( 2 / ( 1 + vrat ) )**pw * vrat**(pw-1)
       
  do irad = 1, nrad
     radius(irad) = rmin * vrat**(float(irad-1)/3.)
     rup(irad) = f1*radius(irad)
     dr(irad) = f2*radius(irad)
  enddo
       
  ! first get the mie coefficients for the condensate
       
  write(miefile,"(A,A,A)")"../Clouds/",trim(cloudname),".mieff"
       
       
  open(unit=10, file= miefile, status='old')
       
       
  read(10,*) idum1, idum2
       
  if (idum1 .ne. nmiewave .or. idum2 .ne. nrad) then
     write(*,*) "Problem with mie coefficients file contents wrong waves or radii  ",trim(miefile)
     stop
  end if
       
  do irad = 1, nrad
     read(10,*) radius_in(irad)
     if (abs(radius(irad) - radius_in(irad)) .gt. &
          0.01*radius(irad)) then
        write(*,*) "Radius grid mismatch in mie file: ",trim(miefile)
        stop
     end if
     do imiewave = 1, nmiewave
        read(10,*) miewavelen(imiewave), qscat(imiewave,irad), &
             qext(imiewave,irad), cos_qscat(imiewave,irad)
     end do
  end do
  close(10)

  
  miewavelen = 1e4 * miewavelen


      ! now for each layer, run through the particle sizes summing up the
    ! cross sections to get optical depth for this cloud from extinction
    ! and scattering
    ! this is all hacked from calc_optics. Credit to Ackerman & Marley
    ext_cloud = 0.d0
    



    if (sizdist .eq. 1) then
       ! log normal

       ! radii supplied in um, convert to cm
       rg = theta(1) * 1e-4
       rsig  = theta(2) * 1e-4
       con = theta(3)
       !  Calculate normalization factor (forces lognormal sum = 1.0)
    
       norm = 0.



       do irad = 1,nrad
          rr = radius(irad)
          arg1 = dr(irad) / ( sqrt(2.*PI)*rr*log(rsig) )
          arg2 = -log( rr/ rg )**2 / ( 2*log(rsig)**2 )
          norm = norm + arg1*exp( arg2 )
       end do
       
       
       ! my constant, CON, is just the column density ndz, in cm-2
       norm = con !/ norm
       
       
       ! now loop over radius and fill up wavelength dependent opacity for
       ! each cloud
       do imiewave = 1, nmiewave
          do irad = 1, 1
             
             rr = radius(irad)
             !arg1 = dr(irad) / ( sqrt(2.*PI)*log(rsig) )
             !arg2 = -log( rr/rg)**2 / ( 2*log(rsig)**2 )
             pir2ndz = norm*PI*rr**2 !*arg1*exp( arg2 )
             
             ext_cloud(imiewave) = ext_cloud(imiewave) + &
                  qext(imiewave,irad)*pir2ndz
          end do ! radius loop
       end do ! miewave loop
       
    else if (sizdist .eq. 2) then
       ! hansen
       ! convert to cm
       a = theta(1) *1e-4
       b = theta(2) *1e-4
       ! already in /cm2
       ndz = theta(3)

       arg1 = ((((2*b) - 1)/b) * log(a*b)) + log(ndz) 
       arg2 = log_gamma((1-(2*b))/b)

       logcon =  (arg1 - arg2) 
       
       do imiewave = 1, nmiewave
          do irad = 1, nrad
             rr = radius(irad)
             drr = dr(irad)
             !write(*,*) (-rr/(a*b)), log(drr)
             arg1 = (-rr/(a*b)) + log(drr)
             !write(*,*) arg1
             arg2 = ((1-3*b)/b) * log(rr)
             arg3 = log(qext(imiewave,irad) * PI * rr**2)
             arg4 = logcon + arg1 + arg2 + arg3
            
             !write(*,*) logcon, arg1, arg2, arg3, arg4 
             ext_cloud(imiewave) = ext_cloud(imiewave) + exp(arg4)
                               
       
          end do
       end do
       
    end if

       
    tau = ext_cloud


  END subroutine cloudtau
