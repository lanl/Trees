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
use infile_variables, only : infuel
use baseline_variables, only : ntspecies,tfuelbins,trhof, &
  grassconstant,litterconstant,lrhof,lfueldepth,lmoist,lsizescale, &
  ldepth,lmoisture,lss,ngrass,gdepth,gmoisture,gss,grho, &
  gmoistoverride
use duet_variables, only : vterminal,StepsPerYear,YearsSinceBurn, &
  Froude,droptime,windprofile,lrhofT,leafdropfreq,decay,grhofT, &
  uavg,vavg,VAR,Umean,Vmean,Uvar,Vvar,fuelSA,lmoistT,gmoistT, &
  lssT,gssT,lafdT,gafdT,compact,moistspec,ssspec,dragco,relhum, &
  periodTotal,litout
use species_variables

implicit none

! Local variables
integer :: i,j,k,ift,ct
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
real,allocatable :: ellarea(:,:),litdecex(:,:,:)
real,dimension(3) :: lmoistsum,lsizesum

! Executable code
call define_duet_variables

trhof = rhof

print*,'Prep work for Duet done'
! For each cell in 3D array of densities, find fall time and x and y 
! displacement for translation of the circle 
do ift=1,ntspecies*tfuelbins
  print*,'Dispersing litter for species',ift
  do i=1,nx
  !print*,'It gets to this point'
    if (MOD(i,int(nx/10)).eq.0) print*,'Placing litter for row',i,'of',nx
    !print*,'It also gets to this point'
    do j=1,ny
    !print*,'ny = ',ny
      do k=1,nz
      !print*,'nz = ',nz
        !print*,'min,max trhof = ',minval(trhof),maxval(trhof)
        !print*,'min, max rhof = ',minval(rhof),maxval(rhof) 
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
                !print*,'inside main loop...'
                if(windprofile.eq.0.or.windprofile.eq.2) then
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
                        !print*,iii,jjj,ellix-xground,ellix+xground,elliy-yground,elliy+yground,ii,jj
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
                !print*,ellarea

                do ii=ellix-xground,ellix+xground
                  if(ii.gt.nx)then
                    ii_real=ii-nx
                  elseif(ii.lt.1)then
                    ii_real=ii+nx
                  else
                    ii_real=ii
                  endif
                  do jj=elliy-yground,elliy+yground
                    if(jj.gt.ny)then
                      jj_real=jj-ny
                    elseif(jj.lt.1)then
                      jj_real=jj+ny
                    else
                      jj_real=jj
                    endif
                    if(ellarea(ii,jj).le.0) ellarea(ii,jj) = 0
                    !print*,ii,jj,ellix-xground,ellix+xground,elliy-yground,elliy+yground
                    !print*,ift,ii_real,jj_real,yt,SHAPE(lrhofT)
                    lrhofT(ift,ii_real,jj_real,yt) = lrhofT(ift,ii_real,jj_real,yt)+ &
                      ellarea(ii,jj)*trhof(ift,i,j,k)/leafdropfreq(ift)
                  !print*,'lrhof ',lrhofT(ift,ii_real,jj_real,yt),ift,ii_real,jj_real,yt
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

if (litout.eq.1) then
  open (11,file='litperyear.dat',form='unformatted',status='unknown')
  write (11) lrhofT
  close (11)


  open  (1,file='beforediff.dat',  form='unformatted',status='unknown')
  write (1) lrhofT
  close (1)

endif

! Diffusion from each cell to the surrounding cells provided the density is less
! than that within the center cell
do j=1,ny
  do i=1,nx
    yt=1
    do yr=1,YearsSinceBurn
      do spy=1,StepsPerYear
        ct=0
        do jj=j-1,j+1
          do ii=i-1,i+1
            if(sum(lrhofT(:,ii,jj,1:yt))-sum(lrhofT(:,i,j,1:yt)).lt.0) then
              ct=ct+1
            endif
          enddo
        enddo
        if(ct.gt.0) then
          do jj=j-3,j+3
            do ii=i-3,i+3
            !if(sum(lrhofT(:,ii,jj,1:yt))-sum(lrhofT(:,i,j,1:yt)).lt.0) then
                do ift=1,ntspecies*tfuelbins
                  lrhofT(ift,ii,jj,yt) = lrhofT(ift,ii,jj,yt) + &
                     0.5/StepsPerYear*1/ct*sum(lrhofT(ift,i,j,1:yt))
                  lrhofT(ift,i,j,yt) = lrhofT(ift,i,j,yt) - &
                     0.5/StepsPerYear*1/ct*sum(lrhofT(ift,i,j,1:yt))
                enddo
            !endif
            enddo
          enddo
        endif
      enddo
    enddo
  enddo
enddo

if (litout.eq.1) then
  open  (1,file='afterdiff.dat',  form='unformatted',status='unknown')
  write (1) lrhofT
  close (1)
  allocate(litdecex(nx,ny,periodTotal))
  do i=1,nx
    do j=1,ny
      yt=1
      do yr=1,YearsSinceBurn
        do spy=1,StepsPerYear
          litdecex(i,j,yt) = lrhofT(1,i,j,1)*exp(-decay(1)*(yt-1))
          yt=yt+1
        enddo
      enddo
    enddo
  enddo
  open  (1,file='litdecex.dat',  form='unformatted',status='unknown')
  write (1) litdecex
  close (1)
endif


! Decay and compact the litter
do ift=1,ntspecies*tfuelbins

  do i=1,nx
    do j=1,ny
      yt=1
      do yr=1,YearsSinceBurn
        do spy=1,StepsPerYear
          lrhofT(ift,i,j,yt) = lrhofT(ift,i,j,yt)*exp(-decay(ift)*(yt-1))
          lafdT(ift,i,j,yt) = lrhofT(ift,i,j,yt)/(dx*dy)*exp(compact(ift)*(yt-1))
          lmoistT(ift,i,j,yt) = lrhofT(ift,i,j,yt)*moistspec(ift)/100*exp(-2.0*yt-1)
          if (lmoistT(ift,i,j,yt).lt.relhum) lmoistT(ift,i,j,yt)=relhum
          if (lrhofT(ift,i,j,yt).gt.0) lssT(ift,i,j,yt) = ssspec(ift) 
        yt=yt+1
        enddo
      enddo
    enddo
  enddo
enddo

if (litout.eq.1) then
  open (12,file='litperyear_dec.dat',form='unformatted',status='unknown')
  write (12) lrhofT
  close (12)
endif


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
      lfueldepth(ift,i,j)=maxval(lafdT(:,i,j,:))
      lmoist(ift,i,j,1)=sum(lmoistT(ift,i,j,:))*sum(lrhofT(ift,i,j,:))/sum(lrhofT(:,i,j,:))
      lsizescale(ift,i,j,1)=sum(lssT(ift,i,j,:))*sum(lrhofT(ift,i,j,:))/sum(lrhofT(:,i,j,:))
    enddo
    do ift=1,ngrass
      rhof(ift+infuel,i,j,1)=grhofT(ift,i,j,YearsSinceBurn*StepsPerYear)
      fueldepth(ift+infuel,i,j,1)=gdepth(ift)
      moist(ift+infuel,i,j,1)=gmoisture(ift)
      sizescale(ift+infuel,i,j,1)=gss(ift)
    enddo
  enddo
enddo

if (gmoistoverride.ne.0) then
  print*,'gmoistoverride = ',gmoistoverride
  print*,'max value of moisture = ',maxval(lmoist(:,:,:,1)+moist(ngrass+infuel:,:,:,1))
  lmoist(:,:,:,1) = gmoistoverride/maxval(lmoist(:,:,:,1)+moist(ngrass+infuel:,:,:,1)) * lmoist(:,:,:,1)
  moist(ngrass+infuel:,:,:,1) = gmoistoverride/maxval(lmoist(:,:,:,1)+moist(ngrass+infuel:,:,:,1)) * moist(ngrass+infuel:,:,:,1)
  print*,'max value of moisture after adjustment = ',maxval(lmoist(:,:,:,1)+moist(ngrass+infuel:,:,:,1))
endif

print*,"Finished DUET litter"
print*,'Litter'
print*,'rhof min/max',minval(lrhof(:,:,:,1)),maxval(lrhof(:,:,:,1))
print*,'afd min/max',minval(lfueldepth(:,:,:)),maxval(lfueldepth(:,:,:))
print*,'mc min/max',minval(lmoist(:,:,:,1)),maxval(lmoist(:,:,:,1))
print*,'ss min/max',minval(lsizescale(:,:,:,1)),maxval(lsizescale(:,:,:,1))
print*,'Grass'
print*,'rhof min/max',minval(rhof(1+infuel:ngrass+infuel,:,:,1)),maxval(rhof(ngrass+infuel:,:,:,1))
print*,'afd min/max',minval(fueldepth(1+infuel:ngrass+infuel,:,:,1)),maxval(fueldepth(ngrass+infuel:,:,:,1))
print*,'mc min/max',minval(moist(1+infuel:ngrass+infuel,:,:,1)),maxval(moist(ngrass+infuel:,:,:,1))
print*,'ss min/max',minval(sizescale(1+infuel:ngrass+infuel,:,:,1)),maxval(sizescale(ngrass+infuel:,:,:,1))

end subroutine Duet
