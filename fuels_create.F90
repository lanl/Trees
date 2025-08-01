!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! fuels_create contains the functions which construct the basic fuel map
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
!-----------------------------------------------------------------------
subroutine fuels_create
!-----------------------------------------------------------------------
! fuels_create is a function which calls the grass and tree 
! fuels_creates and consolidates them to fill the rhof, sizescale, 
! moisture, and fueldepth arrays.
!-----------------------------------------------------------------------
use grid_variables
use infile_variables
use fuels_create_variables
implicit none

! Local variables
integer i,j,k,ift, count
real,external:: zcart

! Executable code
call define_fuels_create_variables 

! Fill tree arrays
if (itrees.ne.0) then
  print*,'Filling Trees fuels_create'
  call tree_fuels_create
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

if(ilitter.eq.2) then
  call TR_RunDUET(nx,ny,nz,ntspecies,ngrass,zheight,trhof,tmoist,tfueldepth,tsizescale,grhof,lrhof, &
  gmoist,lmoist,gsizescale,lsizescale,gfueldepth,lfueldepth)

  print*,'DUET complete.'
endif

! Fill grass arrays
if (igrass.ne.0) then     
  if (ilitter.ne.2) then
    print*,'Filling Grass fuels_create'
    call grass_fuels_create
  endif
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


! Fill litter arrays
if (ilitter.ne.0) then
  if (ilitter.eq.1.or.ilitter.eq.3) then
    if (itrees.gt.0) then
      print*,'Filling Litter fuels_create'
      call litter_fuels_create
    else if (itrees.eq.0) then
      print*,'Warning: itrees=0, no litter placed'
    end if
  !else if (ilitter.eq.2) then
  !  print*,'Filling Litter fuels_create for Duet'
  !  call Duet_Inputs
  endif
  !print*,'ntspecies',ntspecies
  do ift=1,ntspecies
    !print*,'sp: ', ift
    !print*,minval(lrhof(ift,:,:,1)),maxval(lrhof(ift,:,:,1))
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

end subroutine fuels_create

subroutine grass_fuels_create
!-----------------------------------------------------------------------
! grass_fuels_create is a function which computes the characteristics
! of a base grass field and fills rhof, sizescale, moisture, and
! fueldepth arrays.
!-----------------------------------------------------------------------
use grid_variables
use fuels_create_variables
implicit none
logical :: file_exists
! Local variables
integer i,j,k,ift
real target_mass,actual_mass
!real x
real,allocatable:: rhofxy(:,:)

if (igrass.eq.3) then
    ngrass=1
    ! Check if the file exists
    allocate(rhofxy(nx,ny))
    open(1, file="LLM_litter_WG.txt")
     read(1,*) rhofxy
    close(1)
    inquire(file="Saturation.txt", exist=file_exists)
    if (file_exists) then
         ! Read in saturation file to get moisture from ParFlow. AA
      open (1,file='Saturation.txt',form='formatted',status='old')
      print *, 'After opening file??'
      do i=1,10
        do j=1,ny
          do k=1,nx
            read (1,*) satarray(i,j,k)
          enddo
        enddo
      enddo

     print*,'Reading LLM WGlitter'
   endif
endif


! Executable code
do ift=1,ngrass
  do i=1,nx
    do j=1,ny
      if(igrass.eq.3) then 
        grho(ift)=rhofxy(i,j)
        if (file_exists) then
          gmoisture(ift)=satarray(j,i,9)*0.08/0.46
        endif
      endif

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
             
end subroutine grass_fuels_create

subroutine tree_fuels_create 
!-----------------------------------------------------------------------
! tree_fuels_create is a function which computes the characteristics  
! of a forest from the variables designated in 
! define_variables.f. It then fills trhof, tsizescale, tmoisture, 
! and tfueldepth arrays.
!-----------------------------------------------------------------------
use constant_variables
use grid_variables
use fuels_create_variables

implicit none

! Local variables
integer ift,it,i,j,k,ii,jj,kk,iii,jjj,kkk
integer ii_real,jj_real
integer ift_index
!real totarea
real ztest,xloc,yloc
integer xcor,ycor
integer xtop,xbot,ybot,ytop,zbot,ztop
real canopytop,canopybot,canopydiameter,canopymaxh
real atop,abot,test_height,bot_height,top_height
real dbh,barkthick,sstemp,srhoftemp
real tot,add,prt
real target_mass,actual_mass
real,allocatable:: rhoftemp(:)
real,external:: paraboloid,normal 

! Variables for the TreeTracker
integer cellnum, cn
integer cellcount
integer, dimension(100) :: cellid
real, dimension(100) :: cellfuel 
real tfueltot


!real x
! Open tree tracking file:
! open(unit=12,file='TreeTracker.txt')

! Executable Code
!-----Determine the number of trees for each species
print*,'Number of trees of each species:',ntrees

!----- Begin loop which fills arrays with information for each tree
allocate(rhoftemp(tfuelbins))
do i=1,ntspecies
  print*,'Species',i,'with',ntrees(i),'trees'
  do j=1,ntrees(i)
    cellnum = 0
    tfueltot = 0
    if (ntrees(i).gt.9) then
      if (MOD(j,int(ntrees(i)/10)).eq.0) print*,'Placing tree',j,'of',ntrees(i)
    endif
    !----- Place tree location 
    ! Specific tree placement
    if(tlocation(i,j,1).gt.nx*dx.or.tlocation(i,j,1).lt.0) CYCLE
    if(tlocation(i,j,2).gt.ny*dy.or.tlocation(i,j,2).lt.0) CYCLE
    xcor = floor(tlocation(i,j,1)/dx+1)
    ycor = floor(tlocation(i,j,2)/dy+1)


    !----- Determine tree shape characteristics
    ! Shape from tree file
    canopytop = theight(j,i)+zs(xcor,ycor)
    canopybot = tcrownbotheight(j,i)+zs(xcor,ycor)
    canopydiameter = tcrowndiameter(j,i)
    canopymaxh= tcrownmaxheight(j,i)+zs(xcor,ycor)
    
    !----- Translate tree shape to grid
    xbot = floor((tlocation(i,j,1)-canopydiameter/2.)/dx+1)
    xtop = floor((tlocation(i,j,1)+canopydiameter/2.)/dx+1)
    ybot = floor((tlocation(i,j,2)-canopydiameter/2.)/dy+1)
    ytop = floor((tlocation(i,j,2)+canopydiameter/2.)/dy+1)
    
    do k=1,nz-1
      if (canopybot.le.zheight(xcor,ycor,k)) then
        zbot = k
        exit
      endif
    enddo
    do kk=k,nz-1
      if (canopytop.le.zheight(xcor,ycor,kk+1)) then
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
          cellcount = 0
          rhoftemp(:) = 0 ! Density of fuels to be added to current cell of interest
          do iii=1,10
            do jjj=1,10
              do kkk=1,10
                test_height= zheight(ii_real,jj_real,kk)+(2.*kkk-1.)/20. &
                  *(zheight(ii_real,jj_real,kk+1)-zheight(ii_real,jj_real,kk))
                xloc = ((ii-1)+(2.*iii-1.)/20.)*dx
                yloc = ((jj-1)+(2.*jjj-1.)/20.)*dy
                bot_height = paraboloid(abot,xloc,tlocation(i,j,1),yloc,tlocation(i,j,2),canopybot)
                top_height = paraboloid(atop,xloc,tlocation(i,j,1),yloc,tlocation(i,j,2),canopytop)
                if (test_height.ge.bot_height.and.test_height.le.top_height) then 
                  cellcount = cellcount + 1
                  do ift=1,tfuelbins
                    rhoftemp(ift) = rhoftemp(ift)+3./2000.*t2bulkdensity(j,ift,i)* &
                      (test_height-canopybot+4.*(canopytop-canopymaxh)*((xloc-tlocation(i,j,1))**2.+ &
                      (yloc-tlocation(i,j,2))**2.)/canopydiameter**2.)/(canopytop-canopybot) ! Contribution of one subcell to overall bulk density; 
                  enddo
                endif
              enddo
            enddo
          enddo

          !----- Fill in the 3D arrays for a tree
          do ift=1,tfuelbins
            if (cellcount.gt.0) then
              cellnum = cellnum + 1
              cellid(cellnum) = (ii_real+(jj_real*ny)+(kk*ny*nx))
              ! ###### Check the Index here  ######
              cellfuel(cellnum) = rhoftemp(ift)
              tfueltot = tfueltot + cellfuel(cellnum)
            endif
            if (rhoftemp(ift).gt.0) then
              ift_index = (i-1)*ntreefueltypes+ift
              tsizescale(ift_index,ii_real,jj_real,kk) = (trhof(ift_index,ii_real,jj_real,kk)* &
                tsizescale(ift_index,ii_real,jj_real,kk)+rhoftemp(ift)*t2ss(j,ift,i))/ &
                (trhof(ift_index,ii_real,jj_real,kk)+rhoftemp(ift))
              tmoist(ift_index,ii_real,jj_real,kk) = (trhof(ift_index,ii_real,jj_real,kk)* &
                tmoist(ift_index,ii_real,jj_real,kk)+rhoftemp(ift)*t2moisture(j,ift,i))/ &
                (trhof(ift_index,ii_real,jj_real,kk)+rhoftemp(ift))
              trhof(ift_index,ii_real,jj_real,kk) = trhof(ift_index,ii_real,jj_real,kk)+rhoftemp(ift)
            endif
          enddo

        enddo
      enddo
    enddo
    do cn=1, cellnum
       cellfuel(cn) = cellfuel(cn)/tfueltot 
       !print*, cn, cellid(cn), cellfuel(cn), tfueltot
    enddo 
    !Print to TreeTracker.txt file here. AA
    !print*, j, cellnum, cellid(1:cellnum), cellfuel(1:cellnum)
    write(12,*) j, cellnum, cellid(1:cellnum), cellfuel(1:cellnum) 

  enddo
enddo
deallocate(rhoftemp)

!----- Limit Tree Densities for overlapping trees -----
do i=1,nx
  do j=1,ny
    do k=1,zmax
      tot=0.
      do it=1,ntspecies
        ift=trhofmaxindex(it)
        ift_index=(ift-1)*ntreefueltypes
        add=sum(trhof(ift_index+1:ift_index+tfuelbins,i,j,k))
        if(tot+add.le.trhofmax(ift))then
          tot=tot+add
        elseif (add.gt.0)then
          prt=max(0.,trhofmax(ift)-tot)
          trhof(ift_index+1:ift_index+tfuelbins,i,j,k)= &
            trhof(ift_index+1:ift_index+tfuelbins,i,j,k)*(prt/add)
          tot=tot+(prt)
        endif
      enddo
    enddo
  enddo
enddo

!----- SANITY CHECK -----
!----- Tree fuel depth is equal to height of the first cell
do i=1,ntspecies
  do j=1,ntreefueltypes
    ift_index = (i-1)*ntreefueltypes+j
    do ii=1,nx
      do jj=1,ny
        tfueldepth(ift_index,ii,jj) = zheight(ii,jj,2) - zs(ii,jj)
      enddo
    enddo
  enddo
enddo

! Print out the target and actual fuel masses for comparisons sake
target_mass = 0
do i=1,ntspecies
  do j=1,ntrees(i)
    do k=1,tfuelbins
      target_mass = target_mass + t2bulkdensity(j,k,i)*PI*tcrowndiameter(j,i)**2./ &
        8.*(theight(j,i)-tcrownbotheight(j,i))
    enddo
  enddo
enddo
print*,'Trees target fuel mass:',target_mass
actual_mass = 0
do ift=1,ntspecies*ntreefueltypes
  do i=1,nx
    do j=1,ny
      do k=1,zmax
        actual_mass = actual_mass+trhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
      enddo
    enddo
  enddo
enddo
print*,'Trees actual fuel mass:',actual_mass
print*,'Trees error:',100-actual_mass/target_mass*100,'%'

end subroutine tree_fuels_create 

subroutine litter_fuels_create
!-----------------------------------------------------------------------
! litter_fuels_create is a function which computes the characteristics
! of a base litter from trees and fills rhof, sizescale, moisture,
! and fueldepth arrays for litter and subracts mass of the grass 
! arrays for tree shading.
!-----------------------------------------------------------------------
use constant_variables
use grid_variables
use fuels_create_variables
use infile_variables

implicit none

! Local variables
integer i,j,k,ift,ift_grass
real target_mass,actual_mass
real rhocolumn,coverfactor,shadefactor !,rhoftemp
real,allocatable:: rhofxy(:,:)


if (ilitter.eq.3) then
  print*,'Reading LLM tree litter'
  allocate(rhofxy(nx,ny))
  open(1, file="LLM_litter_trees.txt")
    read(1,*) rhofxy
  close(1)

endif

! Executable code
!----- Place litter on ground and remove grass to account for shading
print*,'Placing litter and removing grass to account for shading'
do ift = 1,ntspecies*tfuelbins
  do i=1,nx
    do j=1,ny
      if(ilitter.eq.3) then
        lrho(ift)=rhofxy(i,j)
      endif
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
        if (igrass.eq.1.or.igrass.eq.3) then
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

end subroutine litter_fuels_create 


subroutine tree_baseline 
!-----------------------------------------------------------------
! tree_baseline is a function which computes the characteristics  
! of a forest from the variables designated in 
! define_variables.f. It then fills trhof, tsizescale, tmoisture, 
! and tfueldepth arrays.
!-----------------------------------------------------------------
use constant_variables
use grid_variables
use fuels_create_variables
use tree_tracker_variables

implicit none



integer ift,i,j,k,ii,jj,kk,iii,jjj,kkk
integer ii_real,jj_real
integer ift_index
real totarea
real xtest,ytest
integer xtop,xbot,ybot,ytop,zbot,ztop
real canopytop,canopybot,canopydiameter,canopymaxh
real atop,abot,test_height,bot_height,top_height
integer cellcount
real rhoftemp
real target_mass,actual_mass
real,external:: paraboloid,normal
integer ntrees_sum

integer cellnum, cn
integer, dimension(100) :: cellid
real, dimension(100) :: cellfuel 
real tfueltot
real nsub

nsub = (nx*dx*ny*dy)/(ndatax*ndatay)
!-----Determine the number of trees for each species
if(itrees.eq.1) then
  totarea   = nx*dx*ny*dy
  do i=1,ntspecies
    ntrees(i) = ceiling(totarea*tcanopy(i)/(PI*(tcrowndiametertracker(1,i)/2.)**2.))
  enddo
  print*,'Number of trees of each species:',ntrees
endif
! Open tree tracking file:
open(unit=12,file='TreeTracker.txt')
!----- Begin loop which fills arrays with information for each tree
ntrees_sum = 0
do i=1, ntspecies
  ntrees_sum = ntrees_sum + ntrees(i)
enddo
print*, ntrees
print*, ntrees_sum
do i=1,ntspecies
  if (itrees.eq.1) print*,'Species',i,'with',ntrees(i),'trees'
  print*, 'ADAM ', i, ntspecies, ntrees(i), ntrees(2)
  do j=1, ntrees(i)
    cellnum=0
    tfueltot=0
    if (j.eq.4544) then
      print*, 'ADAM ', j, tlocation(i, j,1), tlocation(i, j,2),tspecies(j)
    endif
    if (MOD(j,1000).eq.0) print*,'Placing tree',j,'of',ntrees(i)
    
    !----- Place tree location
    if (itrees.eq.1.or.itrees.eq.3) then
      ! Randomly place a tree
      call random_number(xtest)
      xtest = xtest*nx*dx
      call random_number(ytest)
      ytest = tdny(1)+ytest*(tdny(2)-tdny(1))*dy
    else
      ! Specific tree placement
      xtest = tlocation(i, j,1)
      ytest = tlocation(i, j,2)
    endif

    ! print*, tspecies
    !----- Determine tree shape characteristics
    if (itrees.eq.1) then
      ! Sample shape from distributions
      canopytop = min(nz*dz,normal(theighttracker(1,i),theighttracker(2,i)))
      canopybot = max(0.,min(canopytop-0.02,normal(tcrownbotheighttracker(1,i),tcrownbotheighttracker(2,i))))
      canopydiameter = max(0.,normal(tcrowndiametertracker(1,i),tcrowndiametertracker(2,i)))
      canopymaxh = min(canopytop-0.01,max(canopybot+0.01,normal(tcrownmaxheighttracker(1,i),tcrownmaxheighttracker(2,i))))
    else
      ! Shape from tree file
      canopytop = theight(j,1)
      canopybot = tcrownbotheight(j,1)
      canopydiameter = tcrowndiameter(j,1)
      canopymaxh= tcrownmaxheight(j,1)
    endif

    !----- Translate tree shape to grid
    xbot = floor((xtest-canopydiameter/2.)/dx+1)
    xtop = floor((xtest+canopydiameter/2.)/dx+1)
    ybot = floor((ytest-canopydiameter/2.)/dy+1)
    ytop = floor((ytest+canopydiameter/2.)/dy+1)
    zbot = 0.
    ztop = 0.
    do k=1,nz-1
      if (SIZE(zheight, DIM=2).le.nint(ytest/dy+1)) then
        zbot = k
        exit
      
      else if (canopybot.le.zheight(nint(xtest/dx+1),nint(ytest/dy+1),k+1)) then
        zbot = k
        exit
      endif
    enddo
    do kk=k,nz-1
      if (SIZE(zheight, DIM=2).le.nint(ytest/dy+1)) then
        zbot = k
        exit
      else if (canopytop.le.zheight(nint(xtest/dx+1),nint(ytest/dy+1),kk+1)) then
        ztop = kk
        if (kk.gt.zmax) zmax=kk
        exit
      endif
    enddo
    
    ! Ellipitical paraboloid used to determine tree contours
    atop = 0.25*(canopydiameter)**2./(canopymaxh-canopytop)
    abot = 0.25*(canopydiameter)**2./(canopymaxh-canopybot)
    if (j.eq.24) then
      print*, "NOLAN HEIGHT STUFF"
      print*, xbot, xtop, zbot, ztop
    endif
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
          cellcount = 0
          do iii=1,10
            do jjj=1,10
              do kkk=1,10
                test_height= zheight(ii_real,jj_real,kk)+(2.*kkk-1.)/20.*(zheight(ii_real,jj_real,kk+1)-zheight(ii_real,jj_real,kk))
                bot_height = paraboloid(abot,((ii-1)+(2.*iii-1.)/20.)*dx,xtest,((jj-1)+(2.*jjj-1.)/20.)*dy,ytest,canopybot)
                top_height = paraboloid(atop,((ii-1)+(2.*iii-1.)/20.)*dx,xtest,((jj-1)+(2.*jjj-1.)/20.)*dy,ytest,canopytop)
                if (test_height.ge.bot_height.and.test_height.le.top_height) then
                  if (j.eq.24) then
                    print*, "NOLAN HEIGHT STUFF", j
                    print*, top_height, bot_height, test_height
                  endif
                  cellcount=cellcount+1
                endif
              enddo
            enddo
          enddo
          if (j.eq.24) then
            print*, "CELL COUNT"
            print*, cellcount
          endif
          !----- Fill in the 3D arrays for a tree
          if (cellcount.gt.0) then
            cellnum = cellnum+1
            do ift=1,tfuelbins
              if (itrees.eq.1) then
                ift_index = (i-1)*tfuelbins+ift
              else
                ift_index = (tspecies(j/nsub +1)-1)*tfuelbins+ift
              endif
              rhoftemp = tbulkdensity(ift,i)*cellcount/1000.
              tsizescale(ift_index,ii_real,jj_real,kk) = (trhof(ift_index,ii_real,jj_real,kk)*tsizescale(ift_index,ii_real,jj_real,kk)+rhoftemp*tss(ift,i))/(trhof(ift_index,ii_real,jj_real,kk)+rhoftemp)
           ! This is the nasty hard code. ~ AA
           ! i==1, LLP, i==2, TurkeyOak
              if(tspecies(j/nsub+1).eq.1)then
               tmoisture(ift,i) = satarray(ii_real,jj_real,4)*1.33/0.55
               !tmoisture(ift,i) = satarray(ii_real,jj_real,4)*1.33/0.55
              elseif(tspecies(j/nsub+1).eq.2)then
               tmoisture(ift,i) = satarray(ii_real,jj_real,6)*1.9/0.5
               !tmoisture(ift,i) = satarray(ii_real,jj_real,5)*2/0.5
              endif
            ! End Hard Code ~ AA
              tmoist(ift_index,ii_real,jj_real,kk) = (trhof(ift_index,ii_real,jj_real,kk)*tmoist(ift_index,ii_real,jj_real,kk)+rhoftemp*tmoisture(ift,i))/(trhof(ift_index,ii_real,jj_real,kk)+rhoftemp)
              trhof(ift_index,ii_real,jj_real,kk) = trhof(ift_index,ii_real,jj_real,kk)+rhoftemp
            enddo
            ! ###### Check the Index here  ###### AA
            ! Domain dimensions are hard coded, and need to be
            ! able to read the domain dimensions
            cellid(cellnum) = (ii_real+(jj_real*ny)+(kk*ny*nx))
            if (j.eq.24) then
              print*, "CELLNUM"
              print*, cellnum
              print*, cellid(cellnum)
            endif
            ! ###### Check the Index here  ######
            cellfuel(cellnum) = rhoftemp
            tfueltot = tfueltot + cellfuel(cellnum)
          endif
        enddo
      enddo
    enddo
    !print*, 'ADAM cellnum ', j, ii_real,jj_real,kk,cellnum
    do cn=1, cellnum
       cellfuel(cn) = cellfuel(cn)/tfueltot 
       !print*, cn, cellid(cn), cellfuel(cn), tfueltot
    enddo 
    !Print to TreeTracker.txt file here. AA
    !print*, j, cellnum, cellid(1:cellnum), cellfuel(1:cellnum)
    if (j.eq.24) then
      print*, "OUTPUT OF ", j
      print*, j, cellnum, cellid(1:cellnum), cellfuel(1:cellnum) 
    endif
    write(12,*) j, cellnum, cellid(1:cellnum), cellfuel(1:cellnum) 
  enddo
  ! if (itrees.ne.1) exit
enddo

close (12) 
end subroutine tree_baseline 

