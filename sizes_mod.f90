module sizes

  implicit none


  ! number of condensate species
  integer, protected :: nclouds

  ! number of gases for mixing
  integer, protected :: ngas

  ! number of patches
  integer,protected :: npatch
  
  ! declares and values size parameters for arrays 
  
  ! number of pressure layers
  integer, protected :: nlayers

  ! max number of layers
  integer,parameter :: maxlayers = 100
  
  ! number of temperatures in line list grid
  !integer, parameter :: nlinetemps = 27
  integer, protected :: nlinetemps
 
  
  ! number of temperatures in CIA tables
  integer,parameter :: nciatemps = 198

  ! number of wavenumbers in lowres CIA tables
  integer,parameter :: ncwave = 1000

  
  ! number of radius bins for mie coeff files
  integer, parameter :: nrad = 60

  ! number of wavelengths in mie coeff files
  integer, parameter :: nmiewave = 196

  ! max wave number
  integer, parameter :: maxwave = 40000

  ! max number of patches
  integer, parameter :: maxpatch = 4
  
  ! number of wavelength/number bins in full range
  integer,protected :: nwave

  save

  
contains
  subroutine initwave(val)
    integer :: val
    nwave = val
  end subroutine initwave

  subroutine initgas(gasval)
    integer :: gasval
    ngas = gasval
  end subroutine initgas

  subroutine initpatch(pval)
    integer :: pval
    if (pval .le. maxpatch) then
       npatch =  pval
    else
       write(*,*) "N patch > maxpatch (", maxpatch,")."
       stop
    end if

  end subroutine initpatch

  subroutine initcloud(pval)
    integer :: pval
    nclouds =  pval
  end subroutine initcloud

  subroutine initlayers(pval)
    integer :: pval
    nlayers =  pval
  end subroutine initlayers

  subroutine inittemps(pval)
    integer :: pval
    nlinetemps =  pval
  end subroutine inittemps
  
end module sizes
