!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! distribution contains functions for number sampling from a distribution
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
real function normal(mu,sigma)
!-----------------------------------------------------------------
! normal is a function which randomly samples a number from a 
! gaussian distribution using the Box-Muller Transform with mu 
! being the mean and sigma the standard deviation
!-----------------------------------------------------------------
use constant_variables

implicit none

! Local Variables
real x1,x2,mu,sigma

! Executable Code
call random_number(x1)
call random_number(x2)
normal = sigma*sqrt(-2.0*log(x1))*cos(2*PI*x2)+mu
 
end function normal
