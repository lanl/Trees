!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! program to populate a grid with an established fuel map
!
! Author: Alexander Josephson (11/19)
! Last Modified: 11/19
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program fuel_maps

      use infile_variables
      use treatment_variables

      print *,'===================================='
      print *,' Running TREES to generate fuel     '
      print *,' files for FIRETEC or QUIC-Fire     '
      print *,'===================================='

      !-----Initialize
      call namelist_input
      call define_constant_variables
      call define_grid_variables
      
      !-----Fuel Read-in
      if(ifuelin.eq.1) call fuel_readin

      !-----Establish baseline
      call baseline

      !-----Perform fuel treatments
      if(itreatment.eq.1) call treatment
      
      !-----Export data to binary files
      call output

      end
