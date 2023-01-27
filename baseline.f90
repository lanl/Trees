!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! baseline contains the functions which construct the basic fuel map
! based off the forest and ground fuel dimensions defined in
! define_variables
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
subroutine baseline
!-----------------------------------------------------------------
! baseline is a function which calls the grass and tree baselines
! and consolidates them to fill the rhof, sizescale, moisture, and
! fueldepth arrays.
!-----------------------------------------------------------------
use grid_variables
use infile_variables
use baseline_variables
implicit none

! Local variables
integer i,j,k,ift
real,external:: zcart

! Executable code
call define_baseline_variables 

! Fill grass arrays
if (igrass.ne.0) then     
  print*,'Filling Grass Baseline'
  call grass_baseline
  do ift=1,ngrass
    do i=1,nx
      do j=1,ny
        do k=1,nz
          rhof(ift+infuel,i,j,k)      = grhof(ift,i,j,k)
          sizescale(ift+infuel,i,j,k) = gsizescale(ift,i,j,k)
          moist(ift+infuel,i,j,k)     = gmoist(ift,i,j,k)
        enddo
        fueldepth(ift+infuel,i,j,1) = gfueldepth(ift,i,j)
      enddo
    enddo
  enddo 
endif

! Fill tree arrays
if (itrees.ne.0) then
  print*,'Filling Trees Baseline'
  call tree_baseline
  do ift=1,ntspecies*ntreefueltypes
    do i=1,nx
      do j=1,ny
        do k=1,nz
          rhof(ift+infuel+ngrass,i,j,k)      = trhof(ift,i,j,k)
          sizescale(ift+infuel+ngrass,i,j,k) = tsizescale(ift,i,j,k)
          moist(ift+infuel+ngrass,i,j,k)     = tmoist(ift,i,j,k)
        enddo
        fueldepth(ift+ngrass+infuel,i,j,1) = tfueldepth(ift,i,j)
      enddo
    enddo
  enddo
endif

! Fill litter arrays
if (ilitter.ne.0) then
  if (ilitter.eq.1) then
    if (itrees.gt.0) then
      print*,'Filling Litter Baseline'
      call litter_baseline
    else if (itrees.eq.0) then
      print*,'Warning: itrees=0, no litter placed'
    end if
  else if (ilitter.eq.2) then
    print*,'Filling Litter Baseline with Duet'
    call Duet
  endif
  do ift=1,ntspecies
    if (sum(trhof(ift,:,:,1)).lt.sum(trhof(ift,:,:,:))*0.01) then
      print*,'Little to no fuel from tree type',ift,'in first cell combining with litter'
      do i=1,nx
        do j=1,ny
          rhof(infuel+ngrass+(ift-1)*ntreefueltypes+1,i,j,1)      = lrhof(ift,i,j,1)
          sizescale(infuel+ngrass+(ift-1)*ntreefueltypes+1,i,j,1) = lsizescale(ift,i,j,1)
          moist(infuel+ngrass+(ift-1)*ntreefueltypes+1,i,j,1)     = lmoist(ift,i,j,1)
          fueldepth(infuel+ngrass+(ift-1)*ntreefueltypes+1,i,j,1) = lfueldepth(ift,i,j)
        enddo
      enddo
    else
      do i=1,nx
        do j=1,ny
          do k=1,nz
            rhof(ift+infuel+ngrass+ntspecies*ntreefueltypes,i,j,k)      = lrhof(ift,i,j,k)
            sizescale(ift+infuel+ngrass+ntspecies*ntreefueltypes,i,j,k) = lsizescale(ift,i,j,k)
            moist(ift+infuel+ngrass+ntspecies*ntreefueltypes,i,j,k)     = lmoist(ift,i,j,k)
          enddo
          fueldepth(ift+infuel+ngrass+ntspecies*ntreefueltypes,i,j,1) = lfueldepth(ift,i,j)
        enddo
      enddo
    endif
  enddo
endif

end subroutine baseline

subroutine grass_baseline
!-----------------------------------------------------------------
! grass_baseline is a function which computes the characteristics
! of a base grass field and fills rhof, sizescale, moisture, and
! fueldepth arrays.
!-----------------------------------------------------------------
use grid_variables
use baseline_variables
implicit none

! Local variables
integer i,j,k,ift
real target_mass,actual_mass
real x

! Executable code
do ift=1,ngrass
  do i=1,nx
    do j=1,ny
      gfueldepth(ift,i,j) = gdepth(ift)
      do k=1,nz-1
        gmoist(ift,i,j,k) = gmoisture(ift)
        gsizescale(ift,i,j,k) = gss(ift)
        if (zheight(i,j,k+1).lt.gdepth(ift)+zs(i,j)) then
          grhof(ift,i,j,k) = grho(ift)
        else
          grhof(ift,i,j,k) = grho(ift)*(gdepth(ift)+zs(i,j)-zheight(i,j,k))/(zheight(i,j,k+1)-zheight(i,j,k))
          if (k.gt.zmax) zmax=k
          exit
        endif
      enddo
    enddo
  enddo
enddo

target_mass = 0
actual_mass = 0
do ift=1,ngrass
  target_mass = target_mass+gdepth(ift)*nx*dx*ny*dy*grho(ift)
  do i=1,nx
    do j=1,ny
      do k=1,zmax
        actual_mass = actual_mass+grhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
      enddo
    enddo
  enddo
enddo
print*,'Grass target fuel mass:',target_mass
print*,'Grass actual fuel mass:',actual_mass
print*,'Grass error:',100-actual_mass/target_mass*100,'%'
             
end subroutine grass_baseline

subroutine tree_baseline 
!-----------------------------------------------------------------
! tree_baseline is a function which computes the characteristics  
! of a forest from the variables designated in 
! define_variables.f. It then fills trhof, tsizescale, tmoisture, 
! and tfueldepth arrays.
!-----------------------------------------------------------------
use constant_variables
use grid_variables
use baseline_variables

implicit none

! Local variables
integer ift,i,j,k,ii,jj,kk,iii,jjj,kkk
integer ii_real,jj_real
integer ift_index
real totarea
real xtest,ytest,ztest,xloc,yloc
integer xcor,ycor
integer xtop,xbot,ybot,ytop,zbot,ztop
real canopytop,canopybot,canopydiameter,canopymaxh
real atop,abot,test_height,bot_height,top_height
real dbh,barkthick,sstemp,srhoftemp
real target_mass,actual_mass
real,allocatable:: rhoftemp(:)
real,external:: paraboloid,normal 
real x

!-----Determine the number of trees for each species
if(itrees.eq.1) then
  totarea = nx*dx*ny*dy
  do i=1,ntspecies
    ntrees(i) = ceiling(totarea*tstemdensity(i)/100.**2.)
  enddo
endif
print*,'Number of trees of each species:',ntrees

!----- Begin loop which fills arrays with information for each tree
allocate(rhoftemp(tfuelbins))
do i=1,ntspecies
  print*,'Species',i,'with',ntrees(i),'trees'
  do j=1,ntrees(i)
    if (ntrees(i).gt.9) then
      if (MOD(j,int(ntrees(i)/10)).eq.0) print*,'Placing tree',j,'of',ntrees(i)
    endif
    !----- Place tree location
    if (itrees.eq.1.or.itrees.eq.3) then
      ! Randomly place a tree
      call random_number(xtest)
      xtest = xtest*nx*dx
      call random_number(ytest)
      ytest = tdny(1)+ytest*(tdny(2)-tdny(1))*dy
    else
      ! Specific tree placement
      xtest = tlocation(i,j,1)
      if(xtest.gt.nx*dx.or.xtest.lt.0) CYCLE
      ytest = tlocation(i,j,2)
      if(ytest.gt.ny*dy.or.ytest.lt.0) CYCLE
    endif
    xcor = floor(xtest/dx+1)
    ycor = floor(ytest/dy+1)

    !----- Determine tree shape characteristics
    if (itrees.eq.1) then
      ! Sample shape from distributions
      canopybot = max(0.,normal(tcrownbotheight(1,i),tcrownbotheight(2,i)))+zs(xcor,ycor)
      canopytop = max(0.002,min(nz*dz,normal(theight(1,i),theight(2,i))))+zs(xcor,ycor)
      canopydiameter = max(0.001,normal(tcrowndiameter(1,i),tcrowndiameter(2,i)))
      canopymaxh = min(canopytop-0.001,max(canopybot+0.001,normal(tcrownmaxheight(1,i),tcrownmaxheight(2,i))+zs(xcor,ycor)))
    else
      ! Shape from tree file
      canopytop = theight(j,i)+zs(xcor,ycor)
      canopybot = tcrownbotheight(j,i)+zs(xcor,ycor)
      canopydiameter = tcrowndiameter(j,i)
      canopymaxh= tcrownmaxheight(j,i)+zs(xcor,ycor)
    endif
    
    !----- Translate tree shape to grid
    xbot = floor((xtest-canopydiameter/2.)/dx+1)
    xtop = floor((xtest+canopydiameter/2.)/dx+1)
    ybot = floor((ytest-canopydiameter/2.)/dy+1)
    ytop = floor((ytest+canopydiameter/2.)/dy+1)
    
    do k=1,nz-1
      if (canopybot.le.zheight(nint(xtest/dx+1),nint(ytest/dy+1),k+1)) then
        zbot = k
        exit
      endif
    enddo
    do kk=k,nz-1
      if (canopytop.le.zheight(nint(xtest/dx+1),nint(ytest/dy+1),kk+1)) then
        ztop = kk
        if (kk.gt.zmax) zmax=kk
        exit
      endif
    enddo
    
    !----- Translate stem and bark to grid
    if(istem.ne.0) then
      ift_index = (i-1)*ntreefueltypes+tfuelbins
      dbh = max(0.001,normal(tdbh(1,i),tdbh(2,i)))
      barkthick = max(0.0001,normal(tbarkthick(1,i),tbarkthick(2,i)))
      do k=1,ztop
        ! Fill stem array
        ztest = min(canopytop,(zheight(xcor,ycor,k+1)+zheight(xcor,ycor,k))/2.)
        sstemp = dbh/2.*(canopytop-ztest)/(canopytop-1.5)
        srhoftemp = trhomicro(i)*PI*sstemp**2./(dy*dx)
        tsizescale(ift_index+1,xcor,ycor,k) = (trhof(ift_index+1,xcor,ycor,k)* &
          tsizescale(ift_index+1,xcor,ycor,k)+srhoftemp*sstemp)/(trhof(ift_index+1,xcor,ycor,k)+srhoftemp)
        tmoist(ift_index+1,xcor,ycor,k) = (trhof(ift_index+1,xcor,ycor,k)* &
          tmoist(ift_index+1,xcor,ycor,k)+srhoftemp*tstemmoist(i))/(trhof(ift_index+1,xcor,ycor,k)+srhoftemp)
        trhof(ift_index+1,xcor,ycor,k) = trhof(ift_index+1,xcor,ycor,k)+srhoftemp

        ! Fill bark array
        srhoftemp = trhomicro(i)*PI*((sstemp+barkthick)**2.-sstemp**2.)/(dy*dx)
        tsizescale(ift_index+2,xcor,ycor,k) = (trhof(ift_index+2,xcor,ycor,k)* &
          tsizescale(ift_index+2,xcor,ycor,k)+srhoftemp*barkthick)/(trhof(ift_index+2,xcor,ycor,k)+srhoftemp)
        tmoist(ift_index+2,xcor,ycor,k) = (trhof(ift_index+2,xcor,ycor,k)* &
          tmoist(ift_index+2,xcor,ycor,k)+srhoftemp*tbarkmoist(i))/(trhof(ift_index+2,xcor,ycor,k)+srhoftemp)
        trhof(ift_index+2,xcor,ycor,k) = trhof(ift_index+2,xcor,ycor,k)+srhoftemp
      enddo
    endif
    
    ! Ellipitical paraboloid used to determine tree canopy contours
    atop = 0.25*canopydiameter**2./(canopymaxh-canopytop)
    abot = 0.25*canopydiameter**2./(canopymaxh-canopybot)
    do ii=xbot,xtop
      if (ii.gt.nx) then
        ii_real = ii-nx 
      else if (ii.lt.1) then 
        ii_real = nx+ii
      else
        ii_real = ii
      endif
      do jj=ybot,ytop
        if (jj.gt.ny) then
          jj_real = jj-ny
        else if (jj.lt.1) then
          jj_real = ny+jj
        else
          jj_real = jj
        endif
        do kk=zbot,ztop
          ! Determine how many of subcells of a cell are within the paraboloid, the fraction of the subcells is equal to the fraction of the cell within the paraboloid
          rhoftemp(:) = 0 ! Density of fuels to be added to current cell of interest
          !print*,'here1'
          do iii=1,10
            do jjj=1,10
              do kkk=1,10
                test_height= zheight(ii_real,jj_real,kk)+(2.*kkk-1.)/20. &
                  *(zheight(ii_real,jj_real,kk+1)-zheight(ii_real,jj_real,kk))
                xloc = ((ii-1)+(2.*iii-1.)/20.)*dx
                yloc = ((jj-1)+(2.*jjj-1.)/20.)*dy
                bot_height = paraboloid(abot,xloc,xtest,yloc,ytest,canopybot)
                top_height = paraboloid(atop,xloc,xtest,yloc,ytest,canopytop)
                if (test_height.ge.bot_height.and.test_height.le.top_height) then 
                  do ift=1,tfuelbins
                    if (itrees.eq.1)then
                      rhoftemp(ift) = rhoftemp(ift)+3./2000.*t1bulkdensity(ift,i)* &
                        (test_height-canopybot+4.*(canopytop-canopymaxh)*((xloc-xtest)**2.+ &
                        (yloc-ytest)**2.)/canopydiameter**2.)/(canopytop-canopybot) ! Contribution of one subcell to overall bulk density
                    else 
                      rhoftemp(ift) = rhoftemp(ift)+3./2000.*t2bulkdensity(j,ift,i)* &
                        (test_height-canopybot+4.*(canopytop-canopymaxh)*((xloc-xtest)**2.+ &
                        (yloc-ytest)**2.)/canopydiameter**2.)/(canopytop-canopybot) ! Contribution of one subcell to overall bulk density; 
                    !note the 4 in the above equations is for making canopydiameter a radius by 
                    !dividing by 2 and then squaring that in the denominator of the denominator
                    endif
                  enddo
                endif
              enddo
            enddo
          enddo

          !----- Fill in the 3D arrays for a tree
          do ift=1,tfuelbins
            if (rhoftemp(ift).gt.0) then
              ift_index = (i-1)*ntreefueltypes+ift
              if (itrees.eq.1) then
                tsizescale(ift_index,ii_real,jj_real,kk) = (trhof(ift_index,ii_real,jj_real,kk)* &
                  tsizescale(ift_index,ii_real,jj_real,kk)+rhoftemp(ift)*t1ss(ift,i))/ &
                  (trhof(ift_index,ii_real,jj_real,kk)+rhoftemp(ift))
                tmoist(ift_index,ii_real,jj_real,kk) = (trhof(ift_index,ii_real,jj_real,kk)* &
                  tmoist(ift_index,ii_real,jj_real,kk)+rhoftemp(ift)*t1moisture(ift,i))/ &
                  (trhof(ift_index,ii_real,jj_real,kk)+rhoftemp(ift))
              else
                tsizescale(ift_index,ii_real,jj_real,kk) = (trhof(ift_index,ii_real,jj_real,kk)* &
                  tsizescale(ift_index,ii_real,jj_real,kk)+rhoftemp(ift)*t2ss(j,ift,i))/ &
                  (trhof(ift_index,ii_real,jj_real,kk)+rhoftemp(ift))
                tmoist(ift_index,ii_real,jj_real,kk) = (trhof(ift_index,ii_real,jj_real,kk)* &
                  tmoist(ift_index,ii_real,jj_real,kk)+rhoftemp(ift)*t2moisture(j,ift,i))/ &
                  (trhof(ift_index,ii_real,jj_real,kk)+rhoftemp(ift))
              endif
            trhof(ift_index,ii_real,jj_real,kk) = trhof(ift_index,ii_real,jj_real,kk)+rhoftemp(ift)
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
deallocate(rhoftemp)

!----- Tree fuel depth is equal to height of the first cell
do i=1,ntspecies
  do j=1,ntreefueltypes
    do ii=tdnx(1),tdnx(2)
      do jj=tdny(1),tdny(2)
        tfueldepth((i-1)*ntreefueltypes+j,ii,jj) = zheight(ii,jj,2)
      enddo
    enddo
  enddo
enddo

! Print out the target and actual fuel masses for comparisons sake
target_mass = 0
if (itrees.eq.1) then
  do i=1,ntspecies
    do j=1,tfuelbins
      target_mass = target_mass + ntrees(i)*t1bulkdensity(j,i)*PI*tcrowndiameter(1,i)**2./ &
        8.*(theight(1,i)-tcrownbotheight(1,i))
    enddo
    if (istem.eq.1) then
      target_mass = target_mass+ntrees(i)*trhomicro(i)*PI*tdbh(1,i)**2.*theight(1,i)/12.
      target_mass = target_mass+ntrees(i)*trhomicro(i)*PI*((tdbh(1,i)+tbarkthick(1,i))**2.-tdbh(1,i)**2.)*theight(1,i)
    endif
  enddo
endif
if (itrees.ne.1) then
  do i=1,ntspecies
    do j=1,ntrees(i)
      do k=1,tfuelbins
        target_mass = target_mass + t2bulkdensity(j,k,i)*PI*tcrowndiameter(j,i)**2./ &
          8.*(theight(j,i)-tcrownbotheight(j,i))
      enddo
    enddo
  enddo
endif
print*,'Trees target fuel mass:',target_mass
actual_mass = 0
do ift=1,ntspecies*ntreefueltypes
  do i=tdnx(1),tdnx(2)
    do j=tdny(1),tdny(2)
      do k=1,zmax
        actual_mass = actual_mass+trhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
      enddo
    enddo
  enddo
enddo
print*,'Trees actual fuel mass:',actual_mass
print*,'Trees error:',100-actual_mass/target_mass*100,'%'

end subroutine tree_baseline 

subroutine litter_baseline
!-----------------------------------------------------------------
! litter_baseline is a function which computes the characteristics
! of a base litter from trees and fills rhof, sizescale, moisture,
! and fueldepth arrays for litter and subracts mass of the grass 
! arrays for tree shading.
!-----------------------------------------------------------------
use constant_variables
use grid_variables
use baseline_variables
use infile_variables

implicit none

! Local variables
integer i,j,k,ift,ift_grass
real target_mass,actual_mass
real rhoftemp,rhocolumn,coverfactor,shadefactor

! Executable code
!----- Place litter on ground and remove grass to account for shading
print*,'Placing litter and removing grass to account for shading'
do ift = 1,ntspecies*tfuelbins
  do i=1,nx
    do j=1,ny
      ! Determine factors for placing litter and removing grass
      rhocolumn = 0
      do k=1,zmax
        rhocolumn = rhocolumn+sum(trhof((ift-1)*tfuelbins+1:ift*tfuelbins,i,j,k))* &
          (zheight(i,j,k+1)-zheight(i,j,k))/theight(ift,1)
      enddo
      if (rhocolumn.gt.0) then
        shadefactor = exp(-grassconstant*rhocolumn/0.6)
        coverfactor = 1.-exp(-litterconstant*rhocolumn/0.6)

        ! Remove grass due to shadefactor
        if (igrass.eq.1) then
          do ift_grass=1,ngrass
            do k=1,nz
              if (zheight(i,j,k).gt.gdepth(ift_grass)+zs(i,j)) then
                exit
              else
                rhof(ift_grass+infuel,i,j,k) = rhof(ift_grass+infuel,i,j,k)*shadefactor
              endif
            enddo
          enddo
        endif

        ! Add litter with dependence to coverfactor
        lfueldepth(ift,i,j) = coverfactor*ldepth(ift)
        do k=1,nz
          if (zheight(i,j,k).gt.lfueldepth(ift,i,j)+zs(i,j)) exit
          if (zheight(i,j,k+1).lt.lfueldepth(ift,i,j)+zs(i,j)) then
            lrhof(ift,i,j,k) = lrho(ift)
          else 
            lrhof(ift,i,j,k) = lrho(ift)*(lfueldepth(ift,i,j)+zs(i,j)-zheight(i,j,k))/(zheight(i,j,k+1)-zheight(i,j,k))
          endif
          lsizescale(ift,i,j,k) = lss(ift)
          lmoist(ift,i,j,k) = lmoisture(ift)
        enddo
      endif
    enddo
  enddo
  print*,'Finished litter for species',ift
enddo

! Print out the target and actual fuel masses for comparisons sake
target_mass = 0
do i=1,ntspecies
  target_mass = target_mass + ntrees(i)*ldepth(i)*PI*tcrowndiameter(1,i)**2.*lrho(i)/4.
enddo
print*,'Litter target fuel mass:',target_mass
actual_mass = 0
do ift=1,ntspecies
  do i=1,nx
    do j=1,ny
      do k=1,zmax
        actual_mass = actual_mass+lrhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
      enddo
    enddo
  enddo
enddo
print*,'Litter actual fuel mass:',actual_mass
print*,'Litter error:',100-actual_mass/target_mass*100,'%'

end subroutine litter_baseline 
