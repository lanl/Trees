!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! metryc contains functions used to define the grid
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real function zcart(aa1,sigma,nz,dz,zsij)
      !-----------------------------------------------------------------
      ! zcart is a function which computes the cartesian vertical 
      ! coordianate
      ! sigma is sigma coordinate
      ! aa1 is the logrithmic stretching coefficient
      ! zs is the vertical height of the ground
      !-----------------------------------------------------------------
      implicit none

      integer nz
      real zb,dz,zsij,f
      real aa1,aa2,aa3    
      real sigma,gdeform
      
      zb = nz*dz
      f  = 0.0
      aa2= f*(1-aa1)/zb
      aa3= (1-aa2*zb-aa1)/zb**2.0
      
      gdeform = aa3*sigma**3.0+aa2*sigma**2.0+aa1*sigma 
      zcart   = gdeform*(zb-zsij)/zb+zsij

      end function zcart
