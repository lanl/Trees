!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! distribution contains functions for number sampling from a distribution
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function normal(mu,sigma)
      !-----------------------------------------------------------------
      ! normal is a function which randomly samples a number from a 
      ! gaussian distribution using the Box-Muller Transform with mu 
      ! being the mean and sigma the standard deviation
      !-----------------------------------------------------------------
      use constant_variables
      
      implicit none

      real x1,x2,mu,sigma
      
      call random_number(x1)
      call random_number(x2)
      normal = sigma*sqrt(-2.0*log(x1))*cos(2*PI*x2)+mu
       
      end function normal
