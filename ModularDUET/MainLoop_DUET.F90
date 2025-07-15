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

! This is a list of subroutines for the DUET program.
! Author:  Jenna Sjunneson McDanold 7/2024


module mainloop
  implicit none
  contains
    !---------------------------------------------------------------------------------------------------
    subroutine disperseLitter
      use DUETio!, only: domain,inarray,outarray,windarray,species
      
      real :: PI = 3.1415926536

      integer :: yr,yt,spy,s,i,ii,iii,j,jj,jjj,k
      integer :: ellix,elliy,xground,yground,ii_real,jj_real
      real :: foliage,zmid,tfall,WalkHorz,WalkSteps,groundradius
      real :: cenx,ceny,ellrotation,magv,magh,inell
      real :: xloc,yloc,untransformedx,untransformedy,dispradius
      real,allocatable :: ellarea(:,:)

      !character*5 :: name = 'MainL' 

      print*,'Subroutine disperseLitter running...'
      print*,'ysb,spy,ns,nx,ny,nz:',domain%ysb,domain%spy,domain%ns,domain%nx,domain%ny,domain%nz

      do yr=1,domain%ysb
        
        do spy=1,domain%spy
          yt=1
          !if(MOD(yt,int(domain%nt/10)).eq.0) print*,'Dispersing litter for year ',yr,' period ',spy
          do s=1,domain%ns
            print*,'Dispersing litter for species',s,'in year',yr
            do i=1,domain%nx
              do j=1,domain%ny
                do k=1,domain%nz-1
                  !print*,'First line = ',inarray%trhof(s,i,j,k),spy,species%step(s)
                  if(inarray%trhof(s,i,j,k).gt.0.and.spy.eq.species%step(s)) then
                    foliage = inarray%trhof(s,i,j,k)
                    !print*,'Foliage = ',foliage
                    zmid = 0.5*(inarray%zheight(i,j,k+1)-inarray%zheight(i,j,k))+inarray%zheight(i,j,k)
                    tfall = zmid/species%vter(s)  ! Time for leaf to hit ground
                    WalkHorz = 0.1*sqrt(species%fuSA(s))/species%frde(s) !! FIXME
                    WalkSteps = ceiling(0.5*tfall)
                    groundradius = sqrt(domain%dx*domain%dy/PI)+WalkSteps*WalkHorz

                      ! Find the center of the ellipse based on the amount of displacement
                      ! from the center of the higrad cell 
                      cenx = i*domain%dx+domain%dx/2.+(windarray%uavg(yt)*tfall)   !rmsum(k,l,2) = xlocs
                      ceny = j*domain%dy+domain%dy/2.+(windarray%vavg(yt)*tfall)   !rmsum(k,l,3) = ylocs   
                      ellix= ceiling(cenx/domain%dx)   ! x index
                      elliy= ceiling(ceny/domain%dy)   ! y index
                    
                      ! Determine the amount of rotation for the ellipse: this is the angle
                      ! between the two displacement vectors (windarray%uavg(yt)*tfall) and (windarray%vavg(yt)*tfall)
                      ellrotation = atan2((windarray%vavg(yt)*tfall),(windarray%uavg(yt)*tfall))                   
                      if(ellrotation.lt.0) ellrotation = ellrotation+2*PI
                    
                      ! Define the horizontal and vertical x and y vectors for stretching the
                      ! ellipse in each direction - horizontal is in perpendicular to the
                      ! displacement vector, and vertical is in line with that vector
                      dispradius = (groundradius*sqrt(zmid**2+(windarray%uavg(yt)*tfall)**2+(windarray%vavg(yt)*tfall)**2))/zmid !r_ground*sqrt(height**2 + |D|**2)/height
                      magv = abs(dispradius+tfall*windarray%uvar(yt))
                      magh = abs(groundradius+tfall*windarray%vvar(yt))
                      xground = ceiling(max(magh,magv)/domain%dx)
                      yground = ceiling(max(magh,magv)/domain%dy)
                      ! Allocate the array ellarea which contains all grid ground cells
                      ! within the area of the ellipse
                      allocate(ellarea(ellix-xground:ellix+xground,elliy-yground:elliy+yground)); ellarea(:,:)=0. !,ELLIPSE(9,(Rx-Lx+1)*(Ry-Ly+1)*100))
                      do ii=ellix-xground,ellix+xground
                        do jj=elliy-yground,elliy+yground
                          ! Calculate proportion of litter deposited in a particular location.
                          inell = 0
                          do iii=1,10
                            xloc = ((ii-1)+(2.*iii-1.)/20.)*domain%dx-cenx
                            do jjj=1,10
                              yloc = ((jj-1)+(2.*jjj-1.)/20.)*domain%dy-ceny
                              untransformedx = xloc*cos(ellrotation)/magh-yloc*sin(ellrotation)/magh
                              untransformedy = xloc*sin(ellrotation)/magv+yloc*cos(ellrotation)/magv
                              if((untransformedx**2.+untransformedy**2.).le.1.) then
                                inell = inell+(1.-(untransformedx**2.+untransformedy**2.)) 
                              endif
                            enddo
                          enddo
                          ellarea(ii,jj) = inell
                        enddo
                      enddo
                      ! Normalize volume in each large grid cell so total is equal to
                      ! volume of 1; find the portion of that volume in each of the grid
                      ! ground cells
                      if(sum(ellarea).lt.1e-8) then
                        if(ellix.lt.1) ellix = ellix + domain%nx
                        if(ellix.gt.domain%nx) ellix = ellix - domain%nx
                        if(elliy.lt.1) elliy = elliy + domain%ny
                        if(elliy.gt.domain%ny) elliy = elliy - domain%ny
                        litter%lrho(s,ellix,elliy,yt) = litter%lrho(s,ellix,elliy,yt)+foliage/species%drop(s)
                      else
                        ellarea = ellarea/sum(ellarea)
                        do ii=ellix-xground,ellix+xground
                          if(ii.gt.domain%nx)then
                            ii_real=ii-domain%nx
                          elseif(ii.lt.1)then
                            ii_real=ii+domain%nx
                          else
                            ii_real=ii
                          endif
                          do jj=elliy-yground,elliy+yground
                            if(jj.gt.domain%ny)then
                              jj_real=jj-domain%ny
                            elseif(jj.lt.1)then
                              jj_real=jj+domain%ny
                            else
                              jj_real=jj
                            endif
                            if(ellarea(ii,jj).le.0) ellarea(ii,jj) = 0
                            litter%lrho(s,ii_real,jj_real,yt) = litter%lrho(s,ii_real,jj_real,yt)+ &
                              ellarea(ii,jj)*foliage/species%drop(s)
                          enddo
                        enddo
                      endif
                      deallocate(ellarea)
                    endif
                enddo
              enddo
            enddo
          enddo
          yt=yt+1
        enddo
      enddo

      print*,'Max and Min of Litter:',maxval(litter%lrho),minval(litter%lrho)
      print*,'Max and Min of trhof:',maxval(inarray%trhof),minval(inarray%trhof)

      print*,'Subroutine disperseLitter complete.'

    end subroutine disperseLitter
    !---------------------------------------------------------------------------------------------------
end module mainloop












