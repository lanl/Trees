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
  use io_variables, only : verbose
  use infile_variables, only : ifuelin
  implicit none
  
  ! Local Variables
  
  ! Executable Code
  print *,'===================================='
  print *,' Running TREES to generate fuel     '
  print *,' files for FIRETEC or QUIC-Fire     '
  print *,'===================================='
  
  !-----Initialize
  call fuellist_input
  call define_grid_variables
  
  !-----Fuel Read-in
  if(ifuelin.eq.1) call grid_readin
  
  !-----Establish fuels_create
  call fuels_create
  
  !-----Export data
  call output_fuel
  if (verbose.eq.1) call output_fuellist

end program fuel_maps
