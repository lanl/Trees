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
  sizescale,moist,fueldepth,trhof,tmoist,tsizescale,tfueldepth
use baseline_variables, only : infuel,ntspecies,tfuelbins, &
  grassconstant,litterconstant,lrhof,lfueldepth,lmoist,lsizescale, &
  ngrass,gdepth,gmoisture,gss,grho, &
  gmoistoverride,fueltotal
use duet_variables, only : vterminal,StepsPerYear,YearsSinceBurn, &
  Froude,droptime,windprofile,lrhofT,leafdropfreq,decay,grhofT, &
  uavg,vavg,VAR,Umean,Vmean,Uvar,Vvar,fuelSA,lmoistT,gmoistT, &
  lssT,gssT,lafdT,gafdT,compact,moistspec,ssspec,dragco,relhum, &
  periodTotal,litout,inputprogram
use species_variables
use FF_variables, only : FFrhof,FFmoist,FFspec,surfrhof,surfdepth,specarray

implicit none

! Local variables
integer :: i,j,k,ift,ct,s,ist,h
integer :: yt,yr,yt3,spy
integer :: ellix,elliy,xground,yground
integer :: ii,jj,ii_real,jj_real
integer :: iii,jjj,fuels
real :: zmid,tfall,totalsum
real :: WalkHorz,WalkSteps,groundradius
real :: xdisp,ydisp,Varx,Vary
real :: cenx,ceny
real :: ellrotation,dispradius
real :: magv,magh,xloc,yloc
real :: untransformedx,untransformedy,inell
real :: rhocolumn,shadefactor,litterFactor
real :: g,foliage
real,allocatable :: ellarea(:,:),litdecex(:,:,:),litsums(:),litprop(:),litprop2(:),specsum(:)
real,dimension(3) :: lmoistsum,lsizesum

! Executable code


if(inputprogram.eq.2) then
  !trhof(1,:,:,:) = FFrhof
  !if(any(isNaN(trhof))) print*,'Reading in NaNs for FF Data rhof'
  !tmoist(1,:,:,:) = FFmoist
  !if(any(isNaN(tmoist))) print*,'Reading in NaNs for FF Data moisture'
  fuels = size(specarray)!+ngrass
  fueltotal = fuels

  deallocate(fuelSA,leafdropfreq,decay, &
  droptime,vterminal,Froude, &
  moistspec,dragco,ssspec,compact, &
  trhof,tmoist,tsizescale,tfueldepth, &
  rhof,moist,sizescale, &
  lrhofT,lmoistT,lssT,lafdT, &
  lrhof,lmoist,lsizescale,lfueldepth, &
  grhofT,gmoistT,gafdT,gssT)

  allocate(fuelSA(fuels),leafdropfreq(fuels),decay(fuels),litprop(fuels), &
  droptime(fuels),vterminal(fuels),Froude(fuels),litsums(fuels), &
  moistspec(fuels),dragco(fuels),ssspec(fuels),compact(fuels),litprop2(fuels), &
  specsum(fuels))

  allocate(rhof(fuels+ngrass,nx,ny,nz)); rhof(:,:,:,:)=0.0
  allocate(sizescale(fuels+ngrass,nx,ny,nz)); sizescale(:,:,:,:)=0.0
  allocate(moist(fuels+ngrass,nx,ny,nz)); moist(:,:,:,:)=0.0
  allocate(fueldepth(fuels+ngrass,nx,ny,nz)); fueldepth(:,:,:,:) = 0.0

  allocate(trhof(fuels+ngrass,nx,ny,nz)); trhof(:,:,:,:)=0.0
  allocate(tsizescale(fuels+ngrass,nx,ny,nz)); tsizescale(:,:,:,:)=0.0
  allocate(tmoist(fuels+ngrass,nx,ny,nz)); tmoist(:,:,:,:)=0.0
  allocate(tfueldepth(fuels+ngrass,nx,ny,nz)); tfueldepth(:,:,:,:)=0.0

  allocate(lrhof(fuels+ngrass,nx,ny,nz)); lrhof(:,:,:,:)=0.0
  allocate(lsizescale(fuels+ngrass,nx,ny,nz)); lsizescale(:,:,:,:)=0.0
  allocate(lmoist(fuels+ngrass,nx,ny,nz)); lmoist(:,:,:,:)=0.0
  allocate(lfueldepth(fuels+ngrass,nx,ny)); lfueldepth(:,:,:)=0.0

  allocate(lrhofT(fuels+ngrass,nx,ny,periodTotal)); lrhofT(:,:,:,:)=0.0
  allocate(lssT(fuels+ngrass,nx,ny,periodTotal)); lssT(:,:,:,:)=0.0
  allocate(lmoistT(fuels+ngrass,nx,ny,periodTotal)); lmoistT(:,:,:,:)=0.0
  allocate(lafdT(fuels+ngrass,nx,ny,periodTotal)); lafdT(:,:,:,:)=0.0

  allocate(grhofT(ngrass,nx,ny,periodTotal)); grhofT(:,:,:,:)=0.0
  allocate(gssT(ngrass,nx,ny,periodTotal)); gssT(:,:,:,:)=0.0
  allocate(gmoistT(ngrass,nx,ny,periodTotal)); gmoistT(:,:,:,:)=0.0
  allocate(gafdT(ngrass,nx,ny,periodTotal)); gafdT(:,:,:,:)=0.0

  do ift = 1,fuels
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (FFspec(i,j,k).eq.specarray(ift)) then
            trhof(ift,i,j,k) = FFrhof(i,j,k)
            tmoist(ift,i,j,k) = FFmoist(i,j,k)
          endif
        enddo
      enddo
    enddo
  enddo

  litsums = 0.0
  litprop = 0.0
  litprop2= 0.0
  totalsum = 0

  do ift=1,fuels
    fuelSA(ift)      =SPECINFO(specarray(ift))%surfarea
    leafdropfreq(ift)=SPECINFO(specarray(ift))%dropperyear
    decay(ift)       =SPECINFO(specarray(ift))%decay
    droptime(ift)    =ceiling(real(SPECINFO(specarray(ift))%stepperyear)/real(12.0/real(StepsPerYear)))
    vterminal(ift)   =SPECINFO(specarray(ift))%vterminal
    Froude(ift)      =SPECINFO(specarray(ift))%froude
    dragco(ift)      =SPECINFO(specarray(ift))%dragco
    moistspec(ift)   =SPECINFO(specarray(ift))%moist
    ssspec(ift)      =SPECINFO(specarray(ift))%sizescale
    compact(ift)     =SPECINFO(specarray(ift))%compact
  enddo

  print*,'leafdropfreq:',leafdropfreq

  allocate(surfrhof(fuels+ngrass,nx,ny),surfdepth(fuels+ngrass,nx,ny))
  print*,'Shape of surfrhof and surfdepth:',shape(surfrhof),shape(surfdepth)
  !fueltotal = fuels
  !allocate(fueldepth(fuels+ngrass,nx,ny,nz)); fueldepth = 0.0
  !print*,'Shape of fueldepth:',shape(fueldepth)

else
  allocate(fueldepth(fueltotal+ngrass,nx,ny,nz)); fueldepth = 0.0
  print*,'Shape of fueldepth:',shape(fueldepth) 
  !fuels = size(specarray)!+ngrass
  !fueltotal = fuels

endif

print*,'fueltotal =',fueltotal
print*,'fuels =',fuels

!allocate(fueldepth(fueltotal+ngrass,nx,ny,nz)); fueldepth = 0.0
!print*,'Shape of fueldepth:',shape(fueldepth)

print*,'max and min of trhof = ',maxval(trhof),minval(trhof)
print*,'Prep work for Duet done'
if(any(isNaN(lrhofT))) print*,'lrhofT is NaN on line 84'
! For each cell in 3D array of densities, find fall time and x and y 
! displacement for translation of the circle 
!if(inputprogram.eq.2) print*,'Dispersing litter'
do ift=1,fueltotal
  print*,'Dispersing litter for species',ift
  if(any(isNaN(lrhofT))) print*,'lrhofT is NaN in line 90'
  do i=1,nx
    do j=1,ny
      do k=1,nz-1
        if(trhof(ift,i,j,k).gt.0) then
          foliage = trhof(ift,i,j,k)
          zmid = 0.5*(zheight(i,j,k+1)-zheight(i,j,k))+zheight(i,j,k)
          !litsums(ift) = litsums(ift) + foliage
          !if (inputprogram.eq.2) then
          !  do s=1,fuels
          !    if (FFSpec(i,j,k).eq.specarray(s)) then
          !      ift = s
          !      litsums(s) = litsums(s) + foliage
          !      exit
          !    endif
          !  enddo
          !  !print*,'ift =',ift,' specarray(ift) = ',specarray(ift),' FFSpec',FFSpec(i,j,k)
          !else 
          !  ift = ist
          !endif
          !print*,'ift=',ift
          !if(any(isNaN(lrhofT))) print*,'lrhofT is NaN in line 156'
          tfall = zmid/vterminal(ift)  ! Time for leaf to hit ground
          !if(isNaN(tfall)) then
          !  print*,'tfall is Nan'
          !  print*,'zmid = ',zmid
          !  print*,'vterminal = ',vterminal(ift)
          !endif
          WalkHorz = 0.1*sqrt(fuelSA(ift))/Froude(ift) !! FIXME
          !if(isNan(WalkHorz)) then
          !  print*,'WalkHorz is Nan'
          !  print*,'fuelSA = ',fuelSA(ift)
          !  print*,'Froude = ',Froude(ift)
          !endif
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
                if(ellrotation.lt.0) ellrotation = ellrotation+2*PI

                ! Define the horizontal and vertical x and y vectors for stretching the
                ! ellipse in each direction - horizontal is in perpendicular to the
                ! displacement vector, and vertical is in line with that vector
                dispradius = (groundradius*sqrt(zmid**2+xdisp**2+ydisp**2))/zmid !r_ground*sqrt(height**2 + |D|**2)/height
                magv = abs(dispradius+tfall*Varx)
                magh = abs(groundradius+tfall*Vary)
                xground = ceiling(max(magh,magv)/dx)
                yground = ceiling(max(magh,magv)/dy)
                !if(any(isNaN(lrhofT))) print*,'lrhofT is NaN in line 186'
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
                        if((untransformedx**2.+untransformedy**2.).le.1.) then
                          inell = inell+(1.-(untransformedx**2.+untransformedy**2.)) 
                        endif
                      enddo
                    enddo
                    ellarea(ii,jj) = inell
                    !if(any(isNaN(ellarea))) then
                    !  print*,'ellarea is NaN'
                    !  print*,'untransformedx',untransformedx
                    !  print*,'untransformedy',untransformedy
                    !endif
                  enddo
                enddo
                ! Normalize volume in each large grid cell so total is equal to
                ! volume of 1; find the portion of that volume in each of the grid
                ! ground cells
                !if(any(isNaN(ellarea))) print*,'ellarea is NaN on line 215'
                !if(sum(ellarea).eq.0) then
                !  print*,'ellarea is 0'
                !  print*,'shape of ellarea = ',shape(ellarea)
                !  print*,'magh, magv = ',magh,magv
                !  print*,'disp and ground radii = ',dispradius,groundradius
                !  print*,''
                !  pause
                !endif
                !if(sum(ellarea).ne.0) ellarea = ellarea/sum(ellarea)
                !if(any(isNaN(lrhofT))) print*,'lrhofT is NaN in line 218'
                !if(any(isNaN(ellarea))) then
                !  print*,'ellarea is NaN on line 219'
                !  print*,'trhof(ift,i,j,k)',trhof(ift,i,j,k)
                !endif
                if(sum(ellarea).eq.0) then
                  if(ellix.lt.1) ellix = ellix + nx
                  if(ellix.gt.nx) ellix = ellix - nx
                  if(elliy.lt.1) elliy = elliy + ny
                  if(elliy.gt.ny) elliy = elliy - ny
                  lrhofT(ift,ellix,elliy,yt) = lrhofT(ift,ellix,elliy,yt)+foliage/leafdropfreq(ift)
                else
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
                      if(jj.gt.ny)then
                        jj_real=jj-ny
                      elseif(jj.lt.1)then
                        jj_real=jj+ny
                      else
                        jj_real=jj
                      endif
                      if(ellarea(ii,jj).le.0) ellarea(ii,jj) = 0
                      lrhofT(ift,ii_real,jj_real,yt) = lrhofT(ift,ii_real,jj_real,yt)+ &
                        ellarea(ii,jj)*foliage/leafdropfreq(ift)
                      !print*,'ift, value',ift,lrhofT(ift,ii_real,jj_real,yt)
                    enddo
                  enddo
                endif
                deallocate(ellarea)
              endif
              yt=yt+1
            enddo
          enddo
        endif
      enddo
    enddo
  if (MOD(i,int(nx/10)).eq.0) then
    print*,'Placed litter for row',i,'of',nx
    !if(any(isNaN(lrhofT))) print*,'The Problem is in LrhofT At The Ellipses'
  endif
enddo
enddo

!print*,'Litsums:',litsums

!print*,'Sum of lrhofT per species:'

if(inputprogram.eq.2) then
!  do h=1,fuels
!    do i=1,nx
!      do j=1,ny
!        if(lrhofT(h,i,j,1).ne.0) then
!          !print*,'species',specarray(h),':',lrhofT(h,i,j,k)
!          specsum(h) = specsum(h) + lrhofT(h,i,j,k)
!          !print*,'Max and min of litter placed so far = ',maxval(lrhofT),minval(lrhofT)
!          !do h=1,fuels
!          !  print*,'Sum of litter for species',specarray(h),':',sum(lrhofT(h,:,:,1))
!          !enddo
!        endif
!      enddo
!    enddo
!  enddo
!  print*,'specsum:',specsum
  do h=1,fuels
    !litprop(h) = litsums(h)/leafdropfreq(h)
    print*,'Sum of litter for species',specarray(h),':',sum(lrhofT(h,:,:,1))
    !print*,'Litter sum / drop prop:',litprop(h)
    !print*,'Max and Min for species',specarray(ift),':',maxval(lrhofT(ift,:,:,1)),minval(lrhofT(ift,:,:,1))
    !print*,'Sum of tree for species',specarray(ift),':',litsums(ift)
  enddo
endif

if (litout.eq.1) then
  open (11,file='litperyear.dat',form='unformatted',status='unknown')
  write (11) lrhofT
  close (11)


  open  (1,file='beforediff.dat',  form='unformatted',status='unknown')
  write (1) lrhofT
  close (1)

endif


!if (inputprogram.eq.2) then
!
!  do ift=1,fuels
!    totalsum = totalsum + sum(lrhofT(ift,:,:,:))
!  enddo
!  do ift=1,fuels
!    litprop(ift) = sum(lrhofT(ift,:,:,:))/totalsum
!    litprop2(ift)= litsums(ift)/sum(trhof(:,:,:,:))
!  enddo
!  print*,'litprop:',litprop
!  print*,'litprop2:',litprop2
!  fueltotal = fuels
!  do ift=1,fuels
!    print*,'Sum of litter for species',specarray(ift),':',sum(lrhofT(ift,:,:,1))
!    !print*,'Max and Min for species',specarray(ift),':',maxval(lrhofT(ift,:,:,1)),minval(lrhofT(ift,:,:,1))
!    !print*,'Sum of tree for species',specarray(ift),':',litsums(ift)
!  enddo
!  fueltotal = fuels
!
!endif


print*,'Beginning diffusion...'
!print*,'Sum before diffusion = ',sum(lrhofT)
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
                do ift=1,fueltotal
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
!if(any(isNaN(lrhofT))) print*,'The Problem is in LrhofT In Diffusion'
print*,'Diffusion complete'
!print*,'Sum after diffusion = ',sum(lrhofT)

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



print*,'Decay and compaction occurring...'
! Decay and compact the litter
do ift=1,fueltotal

  do i=1,nx
    do j=1,ny
      yt=1
      do yr=1,YearsSinceBurn
        do spy=1,StepsPerYear
          lrhofT(ift,i,j,yt) = lrhofT(ift,i,j,yt)*exp(-decay(ift)*(yt-1))
          lafdT(ift,i,j,yt) = min((lrhofT(ift,i,j,yt)/(dx*dy)*exp(compact(ift)*(yt-1))),zheight(i,j,2)-zheight(i,j,1))
          lmoistT(ift,i,j,yt) = lrhofT(ift,i,j,yt)*moistspec(ift)/100*exp(-2.0*yt-1)
          if (lmoistT(ift,i,j,yt).lt.relhum) lmoistT(ift,i,j,yt)=relhum
          if (lrhofT(ift,i,j,yt).gt.0) lssT(ift,i,j,yt) = ssspec(ift) 
        yt=yt+1
        enddo
      enddo
    enddo
  enddo
enddo
print*,'Decay and compaction complete'
if (litout.eq.1) then
  open (12,file='litperyear_dec.dat',form='unformatted',status='unknown')
  write (12) lrhofT
  close (12)
endif
!if(any(isNaN(lrhofT))) print*,'The Problem is in LrhofT in Decay and Compaction'
print*,'Grass growing...'
print*,'grass info = '
print*,'grassconstant = ',grassconstant
! Add grass
do ift=1,ngrass
  !print*,'decay = ',decay(ift)
  do i=1,nx
    do j=1,ny
      ! For every cell that includes grass we will determine the height of the
      ! grass, the amount of shade from the canopy, the the amount of cover from
      ! the litter and previous years of grass, and then build the grass density   
      rhocolumn = 0
      do k=1,zmax
        rhocolumn = rhocolumn+sum(trhof(:,i,j,k))* &
          (zheight(i,j,k+1)-zheight(i,j,k))/zheight(i,j,zmax+1)
          !print*,'zheight(i,j,k+1)-zheight(i,j,k) = ',zheight(i,j,k+1)-zheight(i,j,k)
          !print*,'zheight(i,j,zmax+1) = ',zheight(i,j,zmax+1)
      enddo
      shadefactor = exp(-grassconstant*rhocolumn/0.6)
      yt=1
      do yr=1,YearsSinceBurn
        do spy=1,StepsPerYear
          do yt3=1,yt
            litterFactor=exp(-litterconstant*sum(lrhofT(:,i,j,1:yt3))/0.6)
            call random_number(g)
            if(g.lt.0.75) g = g+0.75
            if(g.gt.1) g = 1
            grhofT(ift,i,j,yt)=grhofT(ift,i,j,yt)+ &
              g*grho(ift)*shadeFactor*litterFactor*exp(-decay(ift)*(yt3-1))
          enddo
          yt=yt+1
        enddo
      enddo
    enddo
  enddo
enddo      
print*,'Grass complete'
print*,'Max, min of grass = ',maxval(grhofT),minval(grhofT)
print*,'Total grass = ',sum(grhofT)
if(any(isNaN(grhofT))) print*,'The Problem is in GrhofT in Grass'

!JENNA FIX ME LATER INPUT PROGRAM EQ 2
!do ift=1,fuels
!  print*,'Sum of litter for species',specarray(ift),':',sum(lrhofT(ift,:,:,1))
!  !print*,'Max and Min for species',specarray(ift),':',maxval(lrhofT(ift,:,:,1)),minval(lrhofT(ift,:,:,1))
!  !print*,'Sum of tree for species',specarray(ift),':',litsums(ift)
!enddo

!if(any(isnan(surfrhof))) print*,'NaNs in surfrhof before any values have been placed'
!if(any(isnan(surfdepth))) print*,'NaNs in surfdepth before any values have been placed'

do i=1,nx
  do j=1,ny
    do ift=1,fueltotal
      lrhof(ift,i,j,1)=sum(lrhofT(ift,i,j,:))
      if(sum(lrhofT(:,i,j,:)).eq.0) then
        lfueldepth(ift,i,j) = 0
      else
        lfueldepth(ift,i,j) = maxval(lafdT(ift,i,j,:))*sum(lrhofT(ift,i,j,:))/sum(lrhofT(:,i,j,:))
      endif
      lmoist(ift,i,j,1)=sum(lmoistT(ift,i,j,:))*sum(lrhofT(ift,i,j,:))/sum(lrhofT(:,i,j,:))
      lsizescale(ift,i,j,1)=sum(lssT(ift,i,j,:))*sum(lrhofT(ift,i,j,:))/sum(lrhofT(:,i,j,:))
    enddo
    !fill grass!
    do ift=1,ngrass
      !if(inputprogram.eq.1) then
        rhof(ift,i,j,1)=grhofT(ift,i,j,YearsSinceBurn*StepsPerYear)
        fueldepth(ift,i,j,1)=gdepth(ift)
        moist(ift,i,j,1)=gmoisture(ift)
        sizescale(ift,i,j,1)=gss(ift)
      !elseif(inputprogram.eq.2) then
        !surfrhof(ift,i,j) = surfrhof(ift,i,j) + grhofT(ift,i,j,YearsSinceBurn*StepsPerYear)
        !if(any(isnan(surfrhof))) print*,'NaNs in grhofT'
        !surfdepth(ift,i,j) = surfdepth(ift,i,j) + gdepth(ift)/(ngrass+fueltotal)
        !if(any(isnan(surfdepth))) print*,'NaNs in gdepth'
      !endif
    enddo
    !fill trees!
    !if(inputprogram.eq.1) then
      do ift=1,fueltotal
        do k=1,nz
          rhof(ift+ngrass,i,j,k)      = trhof(ift,i,j,k)     
          fueldepth(ift+ngrass,i,j,k) = tfueldepth(ift,i,j,k)
          moist(ift+ngrass,i,j,k)     = tmoist(ift,i,j,k)    
          sizescale(ift+ngrass,i,j,k) = tsizescale(ift,i,j,k)
        enddo
      enddo
    !endif
    !fill litter!
    do ift=1,fueltotal
      !if(inputprogram.eq.1) then
        rhof(ift+ngrass+fueltotal,i,j,1)     =lrhof(ift,i,j,1)
        fueldepth(ift+ngrass+fueltotal,i,j,1)=lfueldepth(ift,i,j)
        moist(ift+ngrass+fueltotal,i,j,1)    =lmoist(ift,i,j,1)
        sizescale(ift+ngrass+fueltotal,i,j,1)=lsizescale(ift,i,j,1)
      !elseif(inputprogram.eq.2) then
        !surfrhof(ift+ngrass,i,j) = surfrhof(ift+ngrass,i,j) + lrhof(ift,i,j,1)
        !if(any(isnan(lrhof))) print*,'NaNs in lrhof'
        !surfdepth(ift+ngrass,i,j) = surfdepth(ift+ngrass,i,j) + lfueldepth(ift,i,j)/(ngrass+fueltotal)
        !if(any(isnan(lfueldepth))) print*,'NaNs in lfueldepth'
      !endif
    enddo
  enddo
enddo
!print*,'lrhof max =',maxval(lrhof)
!print*,'grhofT max = ',maxval(grhofT)
surfrhof = rhof(:,:,:,1)
print*,'Shape of surfrhof:',shape(surfrhof)
surfdepth = fueldepth(:,:,:,1)
print*,'Shape of surfdepth:',shape(surfdepth)

if (gmoistoverride.ne.0) then
  print*,'gmoistoverride = ',gmoistoverride
  print*,'max value of moisture = ',maxval(lmoist(:,:,:,1)+moist(ngrass+fueltotal:,:,:,1))
  lmoist(:,:,:,1) = gmoistoverride/maxval(lmoist(:,:,:,1)+moist(ngrass+fueltotal:,:,:,1)) * lmoist(:,:,:,1)
  moist(ngrass+fueltotal:,:,:,1) = gmoistoverride/maxval(lmoist(:,:,:,1)+ &
    moist(ngrass+fueltotal:,:,:,1)) * moist(ngrass+fueltotal:,:,:,1)
  print*,'max value of moisture after adjustment = ',maxval(lmoist(:,:,:,1)+moist(ngrass+fueltotal:,:,:,1))
endif

print*,"Finished DUET litter"
if(inputprogram.eq.1) then
  print*,'Litter' 
  print*,'rhof min/max',minval(lrhof(:,:,:,1)),maxval(lrhof(:,:,:,1))
  if(any(isnan(lrhof))) print*,'STOP.  ERROR.  You have NaNs in Litter Rhof'
  print*,'afd min/max',minval(lfueldepth(:,:,:)),maxval(lfueldepth(:,:,:))
  if(any(isnan(lrhof))) print*,'STOP.  ERROR.  You have NaNs in Litter Fuel Depth'
  print*,'mc min/max',minval(lmoist(:,:,:,1)),maxval(lmoist(:,:,:,1))
  if(any(isnan(lrhof))) print*,'STOP.  ERROR.  You have NaNs in Litter Moisture'
  print*,'ss min/max',minval(lsizescale(:,:,:,1)),maxval(lsizescale(:,:,:,1))
  if(any(isnan(lrhof))) print*,'STOP.  ERROR.  You have NaNs in Litter Sizescale'
  print*,'Grass'
  print*,'rhof min/max',minval(rhof(1:ngrass,:,:,1)),maxval(rhof(1:ngrass,:,:,1))
  if(any(isnan(lrhof))) print*,'STOP.  ERROR.  You have NaNs in Rhof'
  print*,'afd min/max',minval(fueldepth(1:ngrass,:,:,1)),maxval(fueldepth(1:ngrass,:,:,1))
  if(any(isnan(lrhof))) print*,'STOP.  ERROR.  You have NaNs in Fuel Depth'
  print*,'mc min/max',minval(moist(1:ngrass,:,:,1)),maxval(moist(1:ngrass,:,:,1))
  if(any(isnan(lrhof))) print*,'STOP.  ERROR.  You have NaNs in Moisture'
  print*,'ss min/max',minval(sizescale(1:ngrass,:,:,1)),maxval(sizescale(1:ngrass,:,:,1))
  if(any(isnan(lrhof))) print*,'STOP.  ERROR.  You have NaNs in Sizescale'
elseif(inputprogram.eq.2) then
  print*,'surface fuel density min/max',minval(surfrhof),maxval(surfrhof)
  if(any(isnan(lrhof))) print*,'STOP.  ERROR.  You have NaNs in Surface Rhof'
  print*,'surface fuel depth min/max',minval(surfdepth),maxval(surfdepth)
  if(any(isnan(lrhof))) print*,'STOP.  ERROR.  You have NaNs in Surface Depth'
endif
end subroutine Duet
