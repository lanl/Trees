!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine will take a given tree canopy field and produce 
! surface fuels. It is an autonomous program that calculates how the 
! litter will fall from each tree, gather on the ground, and decay and 
! compact over time. Then it calculates where the grass can grow based 
! on the density of litter in that area and how much sunlight the area 
! will get.
!
! Author: Jenna McDanold (2/21)
! Last Modified: 4/22 (AJJ)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine Duet
!-----------------------------------------------------------------
! Using canopy small cell info and grid location, calculate fall 
! time, then x and y displacement based on average wind info. Then
! using this info, find the location on the ground for the 
! displacement and horizontal and vertical stretch to place an 
! elliptical area for dispersing the fuel. Makes an elliptical 
! cylinder with volume 1, and find the area of that ellipse within 
! each grod cell containing the ellipse. This becomes the 
! percentage of the litter within the canopy cell that will be 
! placed within that higrad cell. 
!-----------------------------------------------------------------
use constant_variables, only : PI
use grid_variables, only : nx,ny,nz,dx,dy,zheight,zmax,rhof, &
  sizescale,moist,fueldepth
use baseline_variables, only : ntspecies,tfuelbins,trhof, &
  grassconstant,litterconstant,lrhof,lfueldepth,lmoist,lsizescale, &
  ldepth,lmoisture,lss,ngrass,gdepth,gmoisture,gss,grho
use duet_variables, only : vterminal,StepsPerYear,YearsSinceBurn, &
  Froude,droptime,windprofile,lrhofT,leafdropfreq,decay,grhofT, &
  uavg,vavg,VAR,Umean,Vmean,Uvar,Vvar,fuelSA

implicit none

! Local variables
integer :: i,j,k,ift
integer :: yt,yr,yt3,spy
integer :: ellix,elliy,xground,yground
integer :: ii,jj,ii_real,jj_real
integer :: iii,jjj
real :: zmid,tfall
real :: WalkHorz,WalkSteps,groundradius
real :: xdisp,ydisp,Varx,Vary
real :: cenx,ceny
real :: ellrotation,dispradius
real :: magv,magh,xloc,yloc
real :: untransformedx,untransformedy,inell
real :: rhocolumn,shadefactor,litterFactor
real,allocatable :: ellarea(:,:)

! Executable code
call define_duet_variables

print*,'Prep work for Duet done'

! For each cell in 3D array of densities, find fall time and x and y 
! displacement for translation of the circle 
do ift=1,ntspecies*tfuelbins
  print*,'Dispersing litter for species',ift
  do i=1,nx
    if (MOD(i,int(nx/10)).eq.0) print*,'Placing litter for row',i,'of',nx
    do j=1,ny
      do k=1,nz 
        if(trhof(ift,i,j,k).gt.0) then
          zmid = 0.5*(zheight(i,j,k+1)-zheight(i,j,k))+zheight(i,j,k)

          tfall = zmid/vterminal(ift)  ! Time for leaf to hit ground
          WalkHorz = 0.1*sqrt(fuelSA(ift))/Froude(ift) !! FIXME
          WalkSteps = ceiling(0.5*tfall)
          groundradius = sqrt(dx*dy/PI)+WalkSteps*WalkHorz
    
          !For each year since the last burn...
          yt=1
          do yr=1,YearsSinceBurn
            do spy=1,StepsPerYear
              if(spy.eq.droptime(ift)) then
                if(windprofile.eq.0) then
                  xdisp = uavg(yt)*tfall 
                  ydisp = vavg(yt)*tfall
                  Varx = VAR(yt,1)
                  Vary = VAR(yt,2)
                else  !windprofile=1
                  xdisp = Umean(i,j,k)*tfall  !Distance of leaf displacement (x)
                  ydisp = Vmean(i,j,k)*tfall  !Distance of leaf displacement (y)
                  Varx = Uvar(i,j,k)
                  Vary = Vvar(i,j,k)
                endif

                ! Find the center of the ellipse based on the amount of displacement
                ! from the center of the higrad cell 
                cenx = i*dx+dx/2.+xdisp   !rmsum(k,l,2) = xlocs
                ceny = j*dy+dy/2.+ydisp   !rmsum(k,l,3) = ylocs   
                ellix= ceiling(cenx/dx)   ! x index
                elliy= ceiling(ceny/dy)   ! y index

                ! Determine the amount of rotation for the ellipse: this is the angle
                ! between the two displacement vectors xdisp and ydisp
                ellrotation = atan2(ydisp,xdisp)                   
                if(ellrotation.lt.0) ellrotation = ellrotation+2*pi

                ! Define the horizontal and vertical x and y vectors for stretching the
                ! ellipse in each direction - horizontal is in perpendicular to the
                ! displacement vector, and vertical is in line with that vector
                dispradius = (groundradius*sqrt(zmid**2+xdisp**2+ydisp**2))/zmid !r_ground*sqrt(height**2 + |D|**2)/height
                magv = abs(dispradius+tfall*Varx)
                magh = abs(groundradius+tfall*Vary)
                xground = ceiling(max(magh,magv)/dx)
                yground = ceiling(max(magh,magv)/dy)

                ! Allocate the array ellarea which contains all grid ground cells
                ! within the area of the ellipse
                allocate(ellarea(ellix-xground:ellix+xground,elliy-yground:elliy+yground)); ellarea(:,:)=0. !,ELLIPSE(9,(Rx-Lx+1)*(Ry-Ly+1)*100))
                do ii=ellix-xground,ellix+xground
                  do jj=elliy-yground,elliy+yground
                    ! Calculate proportion of litter deposited in a particular location.
                    inell = 0
                    do iii=1,10
                      xloc = ((ii-1)+(2.*iii-1.)/20.)*dx-cenx
                      do jjj=1,10
                        yloc = ((jj-1)+(2.*jjj-1.)/20.)*dy-ceny
                        untransformedx = xloc*cos(ellrotation)/magh-yloc*sin(ellrotation)/magh
                        untransformedy = xloc*sin(ellrotation)/magv+yloc*cos(ellrotation)/magv
                        if((untransformedx**2.+untransformedy**2.).le.1.) &
                          inell = inell+(1.-(untransformedx**2.+untransformedy**2.)) 
                      enddo
                    enddo
                    ellarea(ii,jj) = inell
                  enddo
                enddo

                ! Normalize volume in each large grid cell so total is equal to
                ! volume of 1; find the portion of that volume in each of the grid
                ! ground cells
                ellarea = ellarea/sum(ellarea)
         
                do ii=ellix-xground,ellix+xground
                  if(ii.gt.nx)then
                    ii_real=ii-nx
                  elseif(ii.lt.1)then
                    ii_real=ii+nx
                  else
                    ii_real=ii
                  endif
                  do jj=elliy-yground,elliy+yground
                    if(jj.gt.nx)then
                      jj_real=jj-nx
                    elseif(jj.lt.1)then
                      jj_real=jj+nx
                    else
                      jj_real=jj
                    endif
                    if(ellarea(ii,jj).le.0) ellarea(ii,jj) = 0
!                    lrhofTemp = ellarea(ii,jj)*trhof(ift,i,j,k)/leafdropfreq(ift)
!                    lmoistT(ift,ii_real,jj_real,yt) = (lrhofT(ift,ii_real,jj_real,yt)* &
!                      lmoistT(ift,ii_real,jj_real,yt)+lrhofTemp*moistsp(ift))/ &
!                      (lrhofT(ift,ii_real,jj_real,yt)+lrhofTemp) ! CONSIDER Variations
                    lrhofT(ift,ii_real,jj_real,yt) = lrhofT(ift,ii_real,jj_real,yt)+ &
                      ellarea(ii,jj)*trhof(ift,i,j,k)/leafdropfreq(ift)
                  enddo
                enddo
                deallocate(ellarea)
              endif
              yt=yt+1
            enddo
          enddo
        endif
      enddo
    enddo
  enddo
enddo

! Decay the litter
do ift=1,ntspecies*tfuelbins
  do i=1,nx
    do j=1,ny
      yt=1
      do yr=1,YearsSinceBurn
        do spy=1,StepsPerYear
          lrhofT(ift,i,j,yt) = lrhofT(ift,i,j,yt)*exp(-decay(ift)*(yt-1))
!           lafdT(i,j,n,yt)=lrhoT(i,j,n,yt)/(dx*dy)*exp(compact(ift)*(yt-1))
        yt=yt+1
        enddo
      enddo
    enddo
  enddo
enddo

! Add grass
do ift=1,ngrass
  do i=1,nx
    do j=1,ny
      ! For every cell that includes grass we will determine the height of the
      ! grass, the amount of shade from the canopy, the the amount of cover from
      ! the litter and previous years of grass, and then build the grass density   
      rhocolumn = 0
      do k=1,zmax
        rhocolumn = rhocolumn+sum(trhof(:,i,j,k))* &
          (zheight(i,j,k+1)-zheight(i,j,k))/zheight(i,j,zmax+1)
      enddo
      shadefactor = exp(-grassconstant*rhocolumn/0.6)
      yt=1
      do yr=1,YearsSinceBurn
        do spy=1,StepsPerYear
          do yt3=1,yt
            litterFactor=exp(-litterconstant*sum(lrhofT(:,i,j,1:yt3))/0.6)
            grhofT(ift,i,j,yt)=grhofT(ift,i,j,yt)+ &
              grho(ift)*shadeFactor*litterFactor*exp(-decay(ift)*(yt3-1))
          enddo
          yt=yt+1
        enddo
      enddo
    enddo
  enddo
enddo      

do i=1,nx
  do j=1,ny
    do ift=1,ntspecies*tfuelbins
      lrhof(ift,i,j,1)=sum(lrhofT(ift,i,j,:))
      lfueldepth(ift,i,j)=ldepth(ift) !max(lafdT(ift,i,j,:))
      lmoist(ift,i,j,1)=lmoisture(ift)
      lsizescale(ift,i,j,1)=lss(ift)
    enddo
    do ift=1,ngrass
      rhof(ift,i,j,1)=grhofT(ift,i,j,YearsSinceBurn*StepsPerYear)
      fueldepth(ift,i,j,1)=gdepth(ift)
      moist(ift,i,j,1)=gmoisture(ift)
      sizescale(ift,i,j,1)=gss(ift)
    enddo
  enddo
enddo

print*,"Finished DUET litter"
print*,'Litter'
print*,'rhof min/max',minval(lrhof(:,:,:,1)),maxval(lrhof(:,:,:,1))
print*,'afd min/max',minval(lfueldepth(:,:,:)),maxval(lfueldepth(:,:,:))
print*,'mc min/max',minval(lmoist(:,:,:,1)),maxval(lmoist(:,:,:,1))
print*,'ss min/mac',minval(lsizescale(:,:,:,1)),maxval(lsizescale(:,:,:,1))
print*,'Grass'
print*,'rhof min/max',minval(rhof(1:ngrass,:,:,1)),maxval(rhof(ngrass:,:,:,1))
print*,'afd min/max',minval(fueldepth(1:ngrass,:,:,1)),maxval(fueldepth(ngrass:,:,:,1))
print*,'mc min/max',minval(moist(1:ngrass,:,:,1)),maxval(moist(ngrass:,:,:,1))
print*,'ss min/max',minval(sizescale(1:ngrass,:,:,1)),maxval(sizescale(ngrass:,:,:,1))

end subroutine Duet
