!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! metryc contains functions used to define the grid
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
real function zcart(aa1,sigma,nz,dz,zsij)
!-----------------------------------------------------------------
! zcart is a function which computes the cartesian vertical 
! coordianate
! sigma is sigma coordinate
! aa1 is the logrithmic stretching coefficient
! zs is the vertical height of the ground
!-----------------------------------------------------------------
implicit none

! Local Variables
integer nz
real zb,dz,zsij,f
real aa1,aa2,aa3    
real sigma,gdeform

! Executable Code
zb = nz*dz
f  = 0.0
aa2= f*(1-aa1)/zb
aa3= (1-aa2*zb-aa1)/zb**2.0

gdeform = aa3*sigma**3.0+aa2*sigma**2.0+aa1*sigma 
zcart   = gdeform*(zb-zsij)/zb+zsij

end function zcart
