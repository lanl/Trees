!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! shapes contains functions to define various tree/treatment shapes
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
real function paraboloid(a,x,x0,y,y0,z0)
!----------------------------------------------------------------
! paraboloid is a function which defines a 3D quadric surface
! a is the sole parameter defined for the function
! x0,y0,z0 are are the point of inflection
! x,y are the tested cartesian coordinates
!----------------------------------------------------------------
implicit none

! Local Variables
real a
real x,y
real x0,y0,z0

! Executable Code
paraboloid = z0+((x-x0)**2.+(y-y0)**2.)/a

end function paraboloid
