!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! program to populate a grid with an established fuel map
!
! Author: Alexander Josephson (11/19)
! Last Modified: 3/22
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! © 2022. Triad National Security, LLC. All rights reserved.  This
! program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by Triad
! National Security, LLC for the U.S.  Department of Energy/National
! Nuclear Security Administration. All rights in the program are
! reserved by Triad National Security, LLC, and the U.S. Department of
! Energy/National Nuclear Security Administration. The Government is
! granted for itself and others acting on its behalf a nonexclusive,
! paid-up, irrevocable worldwide license in this material to reproduce,
! prepare derivative works, distribute copies to the public, perform
! publicly and display publicly, and to permit others to do so.
!-----------------------------------------------------------------
program fuel_maps

use grid_variables
use io_variables
use duet_variables
use species_variables
use baseline_variables
use FF_variables

implicit none

real :: timestart, timefinish
!logical :: there

! Executable Code

! Uncomment the below for FastFuels
print *,'Building input list...'
call FFmakeDuetInputs
print *,'Input list complete.'


!-----Initialize
call cpu_time(timestart)

call namelist_input
call define_constant_variables
call define_grid_variables
call define_baseline_variables

if(inputprogram.eq.1) then
  print *,'===================================='
  print *,' Running DUET to generate fuel     '
  print *,' files for FIRETEC or QUIC-Fire     '
  print *,'===================================='
elseif(inputprogram.eq.2) then
  print *,'===================================='
  print *,' Running DUET to generate fuel     '
  print *,' files for Fastfuels     '
  print *,'===================================='
endif

if(inputprogram.eq.1) then
  call define_species_variables
  call define_duet_variables
  call Duet
elseif(inputprogram.eq.2) then
  call define_ff_variables
  call define_3Dspecies_variables
  call define_3Dduet_variables
  call Duet
endif


!-----Export data to binary files
if(inputprogram.eq.1) then
  if(singlefuel.eq.1)then
    call output_1fuel
  else
    call output_nfuel
  endif
elseif(inputprogram.eq.2) then
  call output_FF
endif

call cpu_time(timefinish)
print*,'Domain size = ',nx,' by ',ny,' by ',nz
if(inputprogram.eq.2) print*,'Average canopy density per area = ',sum(FFrhof)/(nx*dx*ny*dy)
if(inputprogram.eq.2) print*,'Wind Direction = ',winddirection
print*,'Min and Max Winds = ',min(minval(uavg),minval(vavg)),max(maxval(uavg),maxval(vavg))
print*,'Number of timesteps run = ',YearsSinceBurn*StepsPerYear
if (inputprogram.eq.2) then
  print*,'------------------------------'
  print*,'Length of Species array:',size(specarray)
  print*,'Shape of rhof: ',shape(surfrhof)
  print*,'------------------------------'
endif
print*,'Time in seconds = ',timefinish-timestart
print*,'Time in minutes = ',(timefinish-timestart)/60.0

end
