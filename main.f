!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! program to populate a grid with an established fuel map
!
! Author: Alexander Josephson (11/19)
! Last Modified: 11/19
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program fuel_maps

      use grid_variables
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

      end
