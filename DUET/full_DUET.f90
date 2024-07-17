!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine will take a given tree canopy field and produce 
! surface fuels. It is an autonomous program that calculates how the 
! litter will fall from each tree, gather on the ground, and decay and 
! compact over time. Then it calculates where the grass can grow based 
! on the density of litter in that area and how much sunlight the area 
! will get.
!
! Author: Jenna McDanold (2/21)
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
  periodTotal,litout,inputprogram,densitythresh
use species_variables
use FF_variables, only : FFrhof,FFmoist,FFspec,surfrhof,surfdepth,surfmoist,specarray

implicit none

! Local variables
integer :: i,j,k,ift,ct,s,ist,h
integer :: yt,yr,yt3,spy
integer :: ellix,elliy,xground,yground
integer :: ii,jj,ii_real,jj_real
integer :: iii,jjj,fuels,div,nums,iplus,iminus,jplus,jminus
real :: zmid,tfall,totalsum
real :: WalkHorz,WalkSteps,groundradius
real :: xdisp,ydisp,Varx,Vary
real :: cenx,ceny,maxss,premax
real :: ellrotation,dispradius
real :: magv,magh,xloc,yloc
real :: untransformedx,untransformedy,inell
real :: rhocolumn,shadefactor,litterFactor
real :: g,foliage,thresh,a,deltat
real,allocatable :: ellarea(:,:),litdecex(:,:,:),litsums(:),litprop(:),litprop2(:),specsum(:)
real,allocatable :: surfrhof2(:,:,:),surfdepth2(:,:,:),surfmoist2(:,:,:),flattrees(:,:),tmp(:,:,:),TEMP(:,:,:)
real,dimension(3) :: lmoistsum,lsizesum

! Executable code

print*,'ngrass =',ngrass

if(inputprogram.eq.1) specarray = FIA

if(inputprogram.eq.2) then
  !trhof(1,:,:,:) = FFrhof
  !if(any(isNaN(trhof))) print*,'Reading in NaNs for FF Data rhof'
  !tmoist(1,:,:,:) = FFmoist
  !if(any(isNaN(tmoist))) print*,'Reading in NaNs for FF Data moisture'
  fuels = size(specarray)!+ngrass
  !print*,'fuels =',fuels
  !print*,'ngrass =',ngrass
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

  allocate(rhof(2*fuels+ngrass,nx,ny,nz)); rhof(:,:,:,:)=0.0
  allocate(sizescale(2*fuels+ngrass,nx,ny,nz)); sizescale(:,:,:,:)=0.0
  allocate(moist(2*fuels+ngrass,nx,ny,nz)); moist(:,:,:,:)=0.0
  allocate(fueldepth(2*fuels+ngrass,nx,ny,nz)); fueldepth(:,:,:,:) = 0.0

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

  allocate(surfrhof(fuels+ngrass,nx,ny),surfdepth(fuels+ngrass,nx,ny),surfmoist(fuels+ngrass,nx,ny))
  surfrhof = 0.0
  surfdepth = 0.0
  surfmoist = 0.0
  print*,'Shape of surfrhof and surfdepth:',shape(surfrhof),shape(surfdepth)
  !fueltotal = fuels
  !allocate(fueldepth(fuels+ngrass,nx,ny,nz)); fueldepth = 0.0
  !print*,'Shape of fueldepth:',shape(fueldepth)

else
  fuels = fueltotal

  allocate(surfrhof(fuels+ngrass,nx,ny),surfdepth(fuels+ngrass,nx,ny),surfmoist(fuels+ngrass,nx,ny))
  surfrhof = 0.0
  surfdepth = 0.0
  surfmoist = 0.0
  allocate(fueldepth(2*fuels+ngrass,nx,ny,nz)); fueldepth = 0.0
  print*,'Shape of fueldepth:',shape(fueldepth) 
  !fuels = size(specarray)!+ngrass
  !fueltotal = fuels

endif

print*,'fueltotal =',fueltotal
!print*,'fuels =',fuels

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
  !if(any(isNaN(lrhofT))) print*,'lrhofT is NaN in line 90'
  do i=1,nx
    do j=1,ny
      do k=1,nz-1
        if(trhof(ift,i,j,k).gt.0) then
          foliage = trhof(ift,i,j,k)
          zmid = 0.5*(zheight(i,j,k+1)-zheight(i,j,k))+zheight(i,j,k)
          tfall = zmid/vterminal(ift)  ! Time for leaf to hit ground
          WalkHorz = 0.1*sqrt(fuelSA(ift))/Froude(ift) !! FIXME
          WalkSteps = ceiling(0.5*tfall)
          groundradius = sqrt(dx*dy/PI)+WalkSteps*WalkHorz
          yt=1
          do yr=1,YearsSinceBurn
            do spy=1,StepsPerYear
              if(spy.eq.droptime(ift)) then
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
  endif
enddo
enddo


print*,'Sum of lrhofT per species:'
do ift=1,fueltotal
  print*,'Species ',ift,'rhof = ',sum(lrhofT(ift,:,:,:))
enddo
print*,'Sum of litter for all years and species: ',sum(lrhofT)

if(inputprogram.eq.2) then
  do h=1,fuels
    print*,'Sum of litter for species',specarray(h),':',sum(lrhofT(h,:,:,:))
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

print*,'Sum of litter before decay: ',sum(lrhofT)


print*,'Decay and compaction occurring...'
! Decay and compact the litter
do ift=1,fueltotal

  do i=1,nx
    do j=1,ny
      yt=1
      do yr=1,YearsSinceBurn
        do spy=1,StepsPerYear
          !lrhofT(ift,i,j,yt) = lrhofT(ift,i,j,yt)*exp(-decay(ift)*(yt-1))
          !lafdT(ift,i,j,yt) = min((lrhofT(ift,i,j,yt)/(dx*dy)*exp(compact(ift)*(yt-1))),zheight(i,j,2)-zheight(i,j,1))
          !lmoistT(ift,i,j,yt) = lrhofT(ift,i,j,yt)*moistspec(ift)/100*exp(-2.0*yt-1)
          lrhofT(ift,i,j,yt) = lrhofT(ift,i,j,yt)*exp(-decay(ift)*(yt))
          lafdT(ift,i,j,yt) = max((lrhofT(ift,i,j,yt)/(dx*dy) &
            *exp(compact(ift)*(yt))),zheight(i,j,2)-zheight(i,j,1))
          lmoistT(ift,i,j,yt) = lrhofT(ift,i,j,yt)*moistspec(ift)/100*exp(-2.0*yt)
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
print*,'Litter max and min:',maxval(lrhofT),minval(lrhofT)
!if(any(isNaN(lrhofT))) print*,'The Problem is in LrhofT in Decay and Compaction'
print*,'Sum of litter after decay: ',sum(lrhofT)

print*,'Sum of litter before diffusion: ',sum(lrhofT)
print*,'Minimum possible thresh value:',sum(lrhofT)/(nx*ny)
print*,'Density threshold is:',densitythresh

if((sum(lrhofT)/(nx*ny)).ge.densitythresh) then
  if(inputprogram.eq.1) then
    print*,'DANGER WILL ROBINSON!'
    print*,'DENSITYTHRESH IS TOO LOW!!!  MUST BE GREATER THAN ',(1.1)*sum(lrhofT)/(nx*ny)
    print*,'If using TREES, adjust in Fuellist.'
    print*,'STOPPING PROGRAM'
    stop
  else
    densitythresh = sum(lrhofT)/(nx*ny)*2
  endif
endif

!thresh = 0.5
div = 5

a = (dx**2)/4
deltat = (dy**2)/4

allocate(tmp(fuels,nx,ny),TEMP(fuels,nx,ny))
tmp = 0.0
TEMP = 0.0

print*,'Beginning diffusion...'
!print*,'Sum before diffusion = ',sum(lrhofT)
! Diffusion from each cell to the surrounding cells provided the density is less
! than that within the center cell
do j=1,ny
  do i=1,nx
    do ift=1,fuels
      tmp(ift,i,j) = tmp(ift,i,j) + sum(lrhofT(ift,i,j,:))
    enddo
  enddo
enddo

print*,'Sum before diffusion:',sum(tmp)
print*,'Maxval of tmp:',maxval(tmp)

premax = maxval(tmp)

if(maxval(tmp).gt.densitythresh) then
  do while (maxval(tmp).gt.densitythresh)
    print*,'Max of tmp',maxval(tmp)
    if(maxval(tmp).gt.premax) then
      print*,'Diffusion failed.  Continuing without diffusing.'
      TEMP = 0.0
      do j=1,ny
        do i=1,nx
          do ift=1,fuels
            TEMP(ift,i,j) = TEMP(ift,i,j) + sum(lrhofT(ift,i,j,:))
          enddo
        enddo
      enddo      
      continue
    endif
    if(any(isNaN(tmp))) then
      print*,'We have a problem in the diffusion...' 
      stop
    endif
    do j=1,ny
      do i=1,nx
        do ift=1,fuels
          do yt=1,StepsPerYear*YearsSinceBurn
            iplus = i+1
            jplus = j+1
            iminus = i-1
            jminus = j-1

            if (i.eq.1) iminus = nx
            if (i.eq.nx) iplus = 1
            if (j.eq.1) jminus = ny
            if (j.eq.ny) jplus = 1

                ! 2D diffusion equation
            TEMP(ift,i,j) = a*deltat*((tmp(ift,iplus,j)-2*tmp(ift,i,j)+tmp(ift,iminus,j))/dx**2 &
              + (tmp(ift,i,jplus)-2*tmp(ift,i,j)+tmp(ift,i,jminus))/dy**2) &
              + tmp(ift,i,j)
            if((isNaN(TEMP(ift,i,j)))) then
              print*,'The problem is here: ift,i,j',ift,i,j 
              stop
            endif
          enddo
        enddo
      enddo
    enddo!if(any(isNaN(lrhofT))) print*,'The Problem is in LrhofT In Diffusion'
    !do j=1,ny
    !  do i=1,nx
    !    do ift=1,fuels
    !      tmp(ift,i,j) = tmp(ift,i,j) + (TEMP(ift,i,j))
    !    enddo
    !  enddo
    !enddo
    tmp = TEMP
  enddo
else
  TEMP = tmp
endif
!lrhofT(:,:,:,1) = TEMP
print*,'Diffusion complete'
print*,'Sum after diffusion = ',sum(TEMP)
print*,'Sum of lrhofT after diffusion = ',sum(lrhofT)


print*,'Grass growing...'
!print*,'grass info = '
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
      !yt=1
      !do yr=1,YearsSinceBurn
      !  do spy=1,StepsPerYear
      !    do yt3=1,yt
      litterFactor=exp(-litterconstant*sum(TEMP(:,i,j))/0.6)
      !litterFactor=exp(-litterconstant*sum(lrhofT(:,i,j,1:yt3))/0.6)
      !print*,'Litterfactor',litterFactor
      !print*,''
      call random_number(g)
      if(g.lt.0.75) g = g+0.75
      if(g.gt.1) g = 1
      grhofT(ift,i,j,YearsSinceBurn*StepsPerYear)=g*grho(ift)*shadeFactor*litterFactor
            !grhofT(ift,i,j,yt)=grhofT(ift,i,j,yt)+ &
            !  g*grho(ift)*shadeFactor*litterFactor*exp(-decay(ift)*(yt3-1))
      !    enddo
      !    yt=yt+1
      !  enddo
      !enddo
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

lrhof(:,:,:,1)=TEMP
maxss = 100
print*,'Min value of sizescale set to 100'
do i=1,nx
  do j=1,ny
    do ift=1,fueltotal
      !lrhof(ift,i,j,1)=TEMP(ift,i,j)!sum(lrhofT(ift,i,j,:))
      if(sum(lrhof(:,i,j,1)).eq.0) then
        lfueldepth(ift,i,j) = 0
        lmoist(ift,i,j,1) = 0
        lsizescale(ift,i,j,1) = 0
      else
        lfueldepth(ift,i,j)   = maxval(lafdT(ift,i,j,:))*(lrhof(ift,i,j,1))/(sum(grhofT(:,i,j,:))+sum(lrhof(:,i,j,1)))
        lmoist(ift,i,j,1)     = maxval(lmoistT(ift,i,j,:))*(lrhof(ift,i,j,1))/(sum(grhofT(:,i,j,:))+sum(lrhof(:,i,j,1)))
        lsizescale(ift,i,j,1) = maxval(lssT(ift,i,j,:))!*(lrhof(ift,i,j,1))/sum(lrhof(:,i,j,1))
        if(lsizescale(ift,i,j,1).lt.maxss.and.lsizescale(ift,i,j,1).ne.0) then
          maxss = lsizescale(ift,i,j,1)
          print*,'new minss:',maxss
        endif

      endif
    enddo
    !fill grass!
    do ift=1,ngrass
      !if(inputprogram.eq.1) then
        rhof(ift,i,j,1)      = grhofT(ift,i,j,YearsSinceBurn*StepsPerYear)
        if(rhof(ift,i,j,1).ne.0) then
          fueldepth(ift,i,j,1) = gdepth(ift)*(grhofT(ift,i,j,YearsSinceBurn*StepsPerYear)/(sum(grhofT(:,i,j,:))+sum(lrhof(:,i,j,1))))
        !fueldepth(ift,i,j,1) = gdepth(ift)*grhofT(ift,i,j,YearsSinceBurn*StepsPerYear)
          moist(ift,i,j,1)     = gmoisture(ift)*grhofT(ift,i,j,YearsSinceBurn*StepsPerYear)
          sizescale(ift,i,j,1) = gss(ift)!*grhofT(ift,i,j,YearsSinceBurn*StepsPerYear)
        endif
        !if(isnan(moist(ift,i,j,1))) print*,'moist is Nan in grass... ift,i,j',ift,i,j
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
          !if(isnan(moist(ift+ngrass,i,j,1))) print*,'moist is Nan in trees... ift,i,j',ift,i,j
        enddo
      enddo
    !endif
    !fill litter!
    do ift=1,fueltotal
      !if(inputprogram.eq.1) then
      !print*,'iteration:',ift+ngrass
        rhof(ift+ngrass+fueltotal,i,j,1)      = lrhof(ift,i,j,1)
        fueldepth(ift+ngrass+fueltotal,i,j,1) = lfueldepth(ift,i,j)
        moist(ift+ngrass+fueltotal,i,j,1)     = lmoist(ift,i,j,1)
        sizescale(ift+ngrass+fueltotal,i,j,1) = lsizescale(ift,i,j,1)
        !if(isnan(moist(ift+ngrass+fueltotal,i,j,1))) print*,'moist is Nan in litter... ift,i,j',ift,i,j
      !elseif(inputprogram.eq.2) then
        !surfrhof(ift+ngrass,i,j) = surfrhof(ift+ngrass,i,j) + lrhof(ift,i,j,1)
        !if(any(isnan(lrhof))) print*,'NaNs in lrhof'
        !surfdepth(ift+ngrass,i,j) = surfdepth(ift+ngrass,i,j) + lfueldepth(ift,i,j)/(ngrass+fueltotal)
        !if(any(isnan(lfueldepth))) print*,'NaNs in lfueldepth'
      !endif
    enddo
  enddo
enddo
print*,'Min value of sizescale for litter is:',maxss
!print*,'moisture sum of grass in writing:',sum(moist(1,:,:,1))
!print*,'moisture sum of tree1 in writing:',sum(moist(2,:,:,1))
!print*,'moisture sum of tree2 in writing:',sum(moist(3,:,:,1))
!print*,'moisture sum of litt1 in writing:',sum(moist(4,:,:,1))
!print*,'moisture sum of litt2 in writing:',sum(moist(5,:,:,1))
!print*,'moisture sum of litt1 in lrhof:',sum(lmoist(1,:,:,1))
!print*,'moisture sum of litt2 in lrhof:',sum(lmoist(2,:,:,1))
!print*,'sum of rhof2 in lrhof:',sum(lrhof(3,:,:,1))

!print*,'lrhof max =',maxval(lrhof)
!print*,'grhofT max = ',maxval(grhofT)
surfrhof(1,:,:) = rhof(1,:,:,1)
surfdepth(1,:,:) = fueldepth(1,:,:,1)
surfmoist(1,:,:) = moist(1,:,:,1)
print*,'Sum of surface fuel in layer 1 is ',sum(surfrhof(1,:,:))

do s=ngrass+1,ngrass+fuels
  surfrhof(s,:,:) = rhof(fuels+s,:,:,1)
  surfdepth(s,:,:) = fueldepth(fuels+s,:,:,1)
  surfmoist(s,:,:) = moist(fuels+s,:,:,1)
  print*,'Sum of surface fuel in layer ',s,' is ',sum(surfrhof(s,:,:))
enddo
print*,'Shape of surfrhof:',shape(surfrhof)
print*,'Shape of surfdepth:',shape(surfdepth)

allocate(surfrhof2(2,nx,ny),surfdepth2(2,nx,ny),surfmoist2(2,nx,ny),flattrees(nx,ny))

surfrhof2  = 0.0
surfdepth2 = 0.0
surfmoist2 = 0.0

surfrhof2(1,:,:) = surfrhof(1,:,:)
surfdepth2(1,:,:) = surfdepth(1,:,:)
surfmoist2(1,:,:) = surfmoist(1,:,:)

do s=1,fuels
  surfrhof2(2,:,:) = surfrhof2(2,:,:) + surfrhof(s+ngrass,:,:)
  surfdepth2(2,:,:) = surfdepth2(2,:,:) + surfdepth(s+ngrass,:,:)
  surfmoist2(2,:,:) = surfmoist2(2,:,:) + surfmoist(s+ngrass,:,:)
enddo

do j=1,ny
  do i=1,nx
    flattrees(i,j) = sum(trhof(:,i,j,:))
  enddo
enddo

print*,'Sum of each layer of surfrhof2 =',sum(surfrhof2(1,:,:)),sum(surfrhof2(2,:,:))
print*,'Max and Min of surfrhof:',maxval(surfrhof),minval(surfrhof)
print*,'Max and Min of surfdepth:',maxval(surfdepth),minval(surfdepth)
print*,'Max and Min of surfmoist:',maxval(surfmoist),minval(surfmoist)
print*,'Max and Min of rhof:',maxval(rhof(:,:,:,1)),minval(rhof(:,:,:,1))
print*,'Max and Min of depth:',maxval(fueldepth(:,:,:,1)),minval(fueldepth(:,:,:,1))
print*,'Max and Min of moist:',maxval(moist(:,:,:,1)),minval(moist(:,:,:,1))

!print*,'Sum of surfrhof:',sum(surfrhof)
!print*,'Sum of surfdepth:',sum(surfdepth)
!print*,'Sum of surfmoist:',sum(surfmoist)
!print*,'Sum of moist:',sum(moist(:,:,:,1))


if (gmoistoverride.ne.0) then
  print*,'gmoistoverride = ',gmoistoverride
  print*,'max value of moisture = ',maxval(lmoist(:,:,:,1)+moist(ngrass+fueltotal:,:,:,1))
  lmoist(:,:,:,1) = gmoistoverride/maxval(lmoist(:,:,:,1)+moist(ngrass+fueltotal:,:,:,1)) * lmoist(:,:,:,1)
  moist(ngrass+fueltotal:,:,:,1) = gmoistoverride/maxval(lmoist(:,:,:,1)+ &
    moist(ngrass+fueltotal:,:,:,1)) * moist(ngrass+fueltotal:,:,:,1)
  print*,'max value of moisture after adjustment = ',maxval(lmoist(:,:,:,1)+moist(ngrass+fueltotal:,:,:,1))
endif

!print*,'1'
open (12,file='surface_rhof_layered.dat',form='unformatted',status='unknown')
write (12) surfrhof
close (12)
!print*,'2'
open (12,file='surface_depth_layered.dat',form='unformatted',status='unknown')
write (12) surfdepth
close (12)
!print*,'3'
open (12,file='alltrees.dat',form='unformatted',status='unknown')
write (12) flattrees
close (12)

open (12,file='canopy.dat',form='unformatted',status='unknown')
write (12) trhof
close (12)

open (12,file='surface_moist_layered.dat',form='unformatted',status='unknown')
write (12) surfmoist
close (12)

open (12,file='surface_rhof.dat',form='unformatted',status='unknown')
write (12) surfrhof2
close (12)
!print*,'2'
open (12,file='surface_depth.dat',form='unformatted',status='unknown')
write (12) surfdepth2
close (12)
!print*,'3'
open (12,file='surface_moist.dat',form='unformatted',status='unknown')
write (12) surfmoist2
close (12)

!print*,'4'
!print*,'specarray = ',specarray
open (1111,file='surface_species.dat',form='formatted',status='unknown')
  do s=1,size(specarray)
    !print*,'Verifying specarray:',specarray(s)
    write (1111,'(I16)') specarray(s)
  enddo
close(1111)

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
