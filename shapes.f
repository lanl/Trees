!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! shapes contains functions to define various tree/treatment shapes
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function paraboloid(a,x,x0,y,y0,z0)
      !----------------------------------------------------------------
      ! paraboloid is a function which defines a 3D quadric surface
      ! a is the sole parameter defined for the function
      ! x0,y0,z0 are are the point of inflection
      ! x,y are the tested cartesian coordinates
      !----------------------------------------------------------------
      implicit none

      real a
      real x,y
      real x0,y0,z0

      paraboloid = z0+((x-x0)**2.+(y-y0)**2.)/a

      end function paraboloid
