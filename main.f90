!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! program to populate a grid with an established fuel map
!
! Author: Alexander Josephson (11/19)
! Last Modified: 4/23
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
use infile_variables
use treatment_variables
use species_variables
use baseline_variables, only : ilitter,command

implicit none

! Local Variables
logical :: DUETexists

! Executable Code
print *,'===================================='
print *,' Running TREES to generate fuel     '
print *,' files for FIRETEC or QUIC-Fire     '
print *,'===================================='

!-----Initialize
call namelist_input
call define_constant_variables
call define_grid_variables

!-----Fuel Read-in
if(ifuelin.eq.1) call grid_readin

!-----Establish baseline
call baseline

!-----Perform fuel treatments
if(itreatment.ne.0) call treatment

!-----Export data to binary files
print*,'Singlefuel',singlefuel
if(singlefuel.eq.1)then
  call output_1fuel
else
  call output_nfuel
endif

!-----Check for and run DUET
if(ilitter.eq.2) then
  inquire(file='DUET', exist=DUETexists)
  if (DUETexists) then
    print *, "Found DUET"
    call system(command)
  else
    print *, "DUET is not included in this release of Trees! Please use ilitter=1 or 0. See README"
  endif
endif

if (verbose.eq.1) call output_1fuellist

end
