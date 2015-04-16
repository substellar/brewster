module disort_params

  use sizes

  implicit none


    logical,parameter :: LAMBER = .true.
  logical,parameter :: PLANK = .true.
  logical,parameter :: USRTAU = .true.
  logical,parameter :: USRANG = .false.
  logical,parameter :: ONLYFL = .true.

  integer, parameter :: IBCND = 0
  integer, parameter :: NPHI = 0
  integer, parameter :: NMOM = 8
  integer, parameter :: NSTR = 8
  integer, parameter :: nlevel = nlayers+1
  integer, parameter :: NUMU = 8 ! same as NSTR
  integer, parameter :: MAXCLY = nlayers + 1
  integer, parameter :: MAXMOM = NMOM
  integer, parameter :: MAXPHI = 3
  integer, parameter :: MAXULV = 1
  integer, parameter :: MAXUMU = 10
  integer, parameter :: NLYR = nlayers
  real, parameter :: FBEAM = 0.0 ! parallel beam at ToA e.g star/BD
  real, parameter :: FISOT = 0.0  ! not sure what this is
  real, parameter :: ALBEDO = 0.01E0
  real, parameter :: TEMIS = 1e-1 ! gives the top boundary a bit emissivity
  real, parameter :: ACCUR  = 0.0
 
  ! now set up where we want the fluxes.
  ! we just want the top of atmosphere to start
  ! can worry about bonus bits for tracking contributions later

  integer,parameter:: ntau = 1
  real,dimension(ntau),parameter :: utau = 0.01
    
  
  save
  
end module disort_params
