subroutine spex(obspec, modspec,Fratio_int)

  implicit none


  !f2py intent(inout) obspec,modspec
  !f2py intent(out) Fratio_int
  
  integer, parameter:: maxwave = 40000  
  double precision,intent(inout) :: obspec(:,:),modspec(:,:)
  double precision,dimension(maxwave),intent(out) :: Fratio_int
  double precision,allocatable,dimension(:):: delta1,delta2,delta
  double precision,allocatable,dimension(:):: gauss,wlmod,Fmod,wlobs
  double precision:: sigma
  integer :: i, nwmod, nobs,nmod

  nobs = size(obspec(1,:))
  nmod = size(modspec(1,:))

    
  allocate(delta1(nobs),delta2(nobs),delta(nobs),wlobs(nobs))
  allocate(gauss(nmod),wlmod(nmod),Fmod(nmod))

  
  Fmod = modspec(2,:)
  wlmod = modspec(1,:)
  wlobs = obspec(1,:)

  delta1 = abs(wlobs - cshift(wlobs,+1))
  delta2 = abs(wlobs - cshift(wlobs,-1))
  delta = 0.5*(delta1 + delta2)

  
  delta(1) = delta(2)
  delta(nobs) = delta(nobs - 1)

  do i = 1, nobs
     sigma = delta(i) * 3.3 / 2.355
     gauss = exp(-(wlmod - wlobs(i))**2 / (2*sigma**2))
     gauss = gauss / sum(gauss)
     Fratio_int(i) = sum(gauss*Fmod)
  end do

  deallocate(delta1,delta2,delta)
  deallocate(gauss,wlmod)

end subroutine spex

subroutine convFWHM(obspec,modspec,fwhm,Fratio_int)

  implicit none


  !f2py intent(inout) modspec,obspec
  !f2py intent(in) fwhm
  !f2py intent(out) Fratio_int
  
  integer, parameter:: maxwave = 40000  
  double precision,intent(inout) :: modspec(:,:),obspec(:,:) 
  double precision,dimension(maxwave),intent(out) :: Fratio_int
  double precision,allocatable,dimension(:):: wlobs,gauss,wlmod,Fmod
  integer :: i, nwmod, nobs,nmod
  double precision :: fwhm, sigma
  
  nobs = size(obspec(1,:))
  nmod = size(modspec(1,:))

  allocate(wlobs(nobs))
  allocate(gauss(nmod),wlmod(nmod),Fmod(nmod))

  Fmod = modspec(2,:)
  wlmod = modspec(1,:)
  wlobs = obspec(1,:)

  do i = 1, nobs
     sigma = fwhm / 2.355

     gauss = exp(-(wlmod-wlobs(i))**2 / (2*sigma**2))
     
     gauss = gauss/ sum(gauss)
     
     Fratio_int(i)= sum(gauss*Fmod)

  end do

  deallocate(wlobs)
  deallocate(gauss,wlmod,Fmod)

end subroutine convFWHM



subroutine convR(obspec,modspec,R,Fratio_int)

  implicit none


  !f2py intent(inout) modspec,obspec
  !f2py intent(in) R
  !f2py intent(out) Fratio_int
  
  integer, parameter:: maxwave = 40000  
  double precision,intent(inout) :: modspec(:,:),obspec(:,:) 
  double precision,dimension(maxwave),intent(out) :: Fratio_int
  double precision,allocatable,dimension(:):: wlobs,gauss,wlmod,Fmod
  integer :: i, nwmod, nobs,nmod
  double precision:: R, sigma
  
  nobs = size(obspec(1,:))
  nmod = size(modspec(1,:))

  allocate(wlobs(nobs))
  allocate(gauss(nmod),wlmod(nmod),Fmod(nmod))
  
  Fmod = modspec(2,:)
  wlmod = modspec(1,:)
  wlobs = obspec(1,:)

  do i = 1, nobs
     !sigma is FWHM / 2.355
     sigma = (wlobs(i) / R) / 2.355
     gauss = exp(-(wlmod-wlobs(i))**2/(2*sigma**2))
     gauss = gauss/ sum(gauss)
     Fratio_int(i)= sum(gauss*Fmod)
  end do

  deallocate(wlobs)
  deallocate(gauss,wlmod,Fmod)



end subroutine convR
