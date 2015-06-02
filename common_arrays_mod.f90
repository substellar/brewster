module common_arrays

  use sizes
  use define_types
  implicit none

  
  ! declares arrays that will be used alot
  
  ! set up patchy atmosphere object, which is an array of patches
   type(a_patch) :: patch(npatch)

   ! the wavenumber array

   double precision, dimension(nwave) :: wavenum, wavelen
   real,dimension(nlinetemps) :: linetemps

  save

  ! set up the pressure layers - these are hard coded to match the line list
  ! pressure intervals to reduce processing. pressures in mbar
contains
  subroutine set_pressure_scale

    ! pressure scale in mbar
    integer:: i,ipatch

    do ipatch= 1, npatch
       patch(ipatch)%atm%index = (/ (i, i = 1, nlayers) /)

       patch(ipatch)%atm%press=&
            [10.**(-2),10.**(-1.875),10.**(-1.75),10.**(-1.625), &
            10.**(-1.5), 10.**(-1.375),10.**(-1.25), 10.**(-1.125), &
            10.**(-1),10.**(-0.875),10.**(-0.75),10.**(-0.625), &
            10.**(-0.5), 10.**(-0.375),10.**(-0.25), 10.**(-0.125), &
            10.**(0.),10.**(0.125),10.**(0.25),10.**(0.375), &
            10.**(0.5), 10.**(0.625),10.**(0.75), 10.**(0.875), &
            10.**(1.),10.**(1.125),10.**(1.25),10.**(1.375), &
            10.**(1.5), 10.**(1.625),10.**(1.75), 10.**(1.875), &
            10.**(2.),10.**(2.125),10.**(2.25),10.**(2.375), &
            10.**(2.5), 10.**(2.625),10.**(2.75), 10.**(2.875), &
            10.**(3.),10.**(3.125),10.**(3.25),10.**(3.375), &
            10.**(3.5), 10.**(3.625),10.**(3.75), 10.**(3.875), &
            10.**(4.),10.**(4.125),10.**(4.25),10.**(4.375), &
            10.**(4.5), 10.**(4.625),10.**(4.75), 10.**(4.875), &
            10.**(5.),10.**(5.125),10.**(5.25),10.**(5.375), &
            10.**(5.5)]

     
       patch(ipatch)%atm%logP = log10(patch(ipatch)%atm%press)
    end do

  end subroutine set_pressure_scale
  
  subroutine set_line_temps


    linetemps =  [75., 100., 110., 120., 130., 140., 150., 160., 170., 180., &
       190., 200., 210., 220., 230., 240., 250., 260., 270., 275., 280., 290.,&
       300., 310., 320., 330., 340., 350., 375., 400., 425., 450., 475., &
       500., 575., 650., 725., 800., 900., 1000., 1100., 1200., 1300., &
       1400., 1500., 1600., 1700., 1800., 1900., 2000., 2150., 2300., &
       2450., 2600., 2800., 3000., 3250., 3500., 3750., 4000.]

    !old line temps
!    [75., 85., 100., 120., 140., 160., 180., 200., 230., 260., &
 !          300., 350., 400.,500., 650., 800., 1000., 1200., 1400., 1600.,&
  !         1800., 2000., 2300., 2600.,3000., 3500., 4000.]



  end subroutine set_line_temps



end module common_arrays
