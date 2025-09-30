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

!-----------------------------------------------------------------------
! fuels_create is a function which calls the grass and tree 
! fuels_creates and consolidates them to fill the rhof, sizescale, 
! moisture, and fueldepth arrays.
!-----------------------------------------------------------------------
subroutine fuels_create
  use grid_variables, only : nx,ny,nz,rhof,sizescale,moist,fueldepth, &
    zheight
  use infile_variables, only : infuel
  use fuel_variables, only : ntspecies,ntreefueltypes,trhof,ilitter, &
    tsizescale,tmoist,tfueldepth,ngrass,grhof,lrhof,gmoist,lmoist, &
    gsizescale,lsizescale,gfueldepth,lfueldepth,itrees,igrass
  implicit none
  
  ! Local variables
  integer i,j,k,ift,it
  real,external:: zcart
  
  ! Executable code
  call define_fuel_variables 
  
  ! Fill tree arrays
  if (itrees.ne.0) then
    print*,'Filling Trees fuels_create'
    call tree_fuels_create
    do ift=1,ntspecies*ntreefueltypes
      it=ift+infuel+ngrass
      do i=1,nx
        do j=1,ny
          do k=1,nz
            rhof(it,i,j,k)      = trhof(ift,i,j,k)
            sizescale(it,i,j,k) = tsizescale(ift,i,j,k)
            moist(it,i,j,k)     = tmoist(ift,i,j,k)
          enddo
          fueldepth(it,i,j,1) = tfueldepth(ift,i,j)
        enddo
      enddo
    enddo
  endif
  
  ! Fill grass arrays
  if (igrass.ne.0) then     
    if (ilitter.ne.2) then
      print*,'Filling Grass fuels_create'
      call grass_fuels_create
    endif
    do ift=1,ngrass
      it=ift+infuel
      do i=1,nx
        do j=1,ny
          do k=1,nz
            rhof(it,i,j,k)      = grhof(ift,i,j,k)
            sizescale(it,i,j,k) = gsizescale(ift,i,j,k)
            moist(it,i,j,k)     = gmoist(ift,i,j,k)
          enddo
          fueldepth(it,i,j,1) = gfueldepth(ift,i,j)
        enddo
      enddo
    enddo 
  endif
  
  
  ! Fill litter arrays
  if (ilitter.ne.0) then
    if (ilitter.eq.1) then
      if (itrees.gt.0) then
        print*,'Filling Litter fuels_create'
        call litter_fuels_create
      else if (itrees.eq.0) then
        print*,'Warning: itrees=0, no litter placed'
      end if
    else if (ilitter.eq.2) then
      call TR_RunDUET(nx,ny,nz,ntspecies,ngrass,zheight,trhof,tmoist, &
        tfueldepth,tsizescale,grhof,lrhof,gmoist,lmoist,gsizescale, &
        lsizescale,gfueldepth,lfueldepth)
      print*,'DUET complete.'
    endif
    do ift=1,ntspecies
      if (sum(trhof(ift,:,:,1)).lt.sum(trhof(ift,:,:,:))*0.01) then
        print*,'Little to no fuel from tree type',ift, &
          'in first cell combining with litter'
        it=infuel+ngrass+(ift-1)*ntreefueltypes+1
        do i=1,nx
          do j=1,ny
            rhof(it,i,j,1)      = lrhof(ift,i,j,1)
            sizescale(it,i,j,1) = lsizescale(ift,i,j,1)
            moist(it,i,j,1)     = lmoist(ift,i,j,1)
            fueldepth(it,i,j,1) = lfueldepth(ift,i,j)
          enddo
        enddo
      else
        it=ift+infuel+ngrass+ntspecies*ntreefueltypes
        do i=1,nx
          do j=1,ny
            do k=1,nz
              rhof(it,i,j,k)      = lrhof(ift,i,j,k)
              sizescale(it,i,j,k) = lsizescale(ift,i,j,k)
              moist(it,i,j,k)     = lmoist(ift,i,j,k)
            enddo
            fueldepth(it,i,j,1) = lfueldepth(ift,i,j)
          enddo
        enddo
      endif
    enddo
  endif

end subroutine fuels_create

!-----------------------------------------------------------------------
! grass_fuels_create is a function which computes the characteristics
! of a base grass field and fills rhof, sizescale, moisture, and
! fueldepth arrays.
!-----------------------------------------------------------------------
subroutine grass_fuels_create
  use grid_variables, only : nx,ny,nz,dx,dy,zheight,zs
  use fuel_variables, only : ngrass,gfueldepth,gdepth,gmoist, &
    gmoisture,gsizescale,gss,grhof,grho
  implicit none
  
  ! Local variables
  integer :: i,j,k,ift
  integer :: zmax
  real :: target_mass,actual_mass
  
  ! Executable code
  zmax=0
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
            grhof(ift,i,j,k) = grho(ift)* &
              (gdepth(ift)+zs(i,j)-zheight(i,j,k))/ &
              (zheight(i,j,k+1)-zheight(i,j,k))
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
          actual_mass = actual_mass+ &
            grhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
        enddo
      enddo
    enddo
  enddo
  print*,'Grass target fuel mass:',target_mass
  print*,'Grass actual fuel mass:',actual_mass
  print*,'Grass error:',100-actual_mass/target_mass*100,'%'
end subroutine grass_fuels_create

!-----------------------------------------------------------------------
! tree_fuels_create is a function which computes the characteristics  
! of a forest from the variables designated in 
! define_variables.f. It then fills trhof, tsizescale, tmoisture, 
! and tfueldepth arrays.
!-----------------------------------------------------------------------
subroutine tree_fuels_create
  use constant_variables, only : PI 
  use grid_variables, only : nx,ny,nz,dx,dy,zs,zheight,zmax
  use fuel_variables, only : tfuelbins,ntspecies,ntrees,tlocation, &
    theight,tcrownbotheight,tcrowndiameter,tcrownmaxheight,tfuelbins, &
    ntreefueltypes,istem,tbarkthick,tsizescale,trhof,tmoist,tfueldepth,&
    t2bulkdensity,tdbh,trhomicro,tstemmoist,tbarkmoist,t2ss,t2moisture,&
    trhofmaxindex,trhofmax
  implicit none
  
  ! Local variables
  integer :: ift,it,i,j,k,ii,jj,kk,iii,jjj,kkk
  integer :: ii_real,jj_real
  integer :: ift_index
  integer :: xcor,ycor
  integer :: xtop,xbot,ybot,ytop,zbot,ztop
  real :: ztest,xloc,yloc
  real :: canopytop,canopybot,canopydiameter,canopymaxh
  real :: atop,abot,test_height,bot_height,top_height
  real :: dbh,barkthick,sstemp,srhoftemp
  real :: tot,add,prt
  real :: target_mass,actual_mass
  real,allocatable :: rhoftemp(:)
  real,external :: paraboloid,normal 
  
  ! Executable Code
  !-----Determine the number of trees for each species
  print*,'Number of trees of each species:',ntrees
  
  !----- Begin loop which fills arrays with information for each tree
  allocate(rhoftemp(tfuelbins))
  do i=1,ntspecies
    print*,'Species',i,'with',ntrees(i),'trees'
    do j=1,ntrees(i)
      if (ntrees(i).gt.9) then
        if (MOD(j,int(ntrees(i)/10)).eq.0) print*, &
          'Placing tree',j,'of',ntrees(i)
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
          ztest = min(canopytop,(zheight(xcor,ycor,k+1)+ &
            zheight(xcor,ycor,k))/2.)
          sstemp = dbh/2.*(canopytop-ztest)/(canopytop-1.5)
          srhoftemp = trhomicro(i)*PI*sstemp**2./(dy*dx)
          tsizescale(ift_index+1,xcor,ycor,k) = &
            (trhof(ift_index+1,xcor,ycor,k)* &
            tsizescale(ift_index+1,xcor,ycor,k)+srhoftemp*sstemp)/ &
            (trhof(ift_index+1,xcor,ycor,k)+srhoftemp)
          tmoist(ift_index+1,xcor,ycor,k) = &
            (trhof(ift_index+1,xcor,ycor,k)* &
            tmoist(ift_index+1,xcor,ycor,k)+srhoftemp*tstemmoist(i))/ &
            (trhof(ift_index+1,xcor,ycor,k)+srhoftemp)
          trhof(ift_index+1,xcor,ycor,k) = &
            trhof(ift_index+1,xcor,ycor,k)+srhoftemp
  
          ! Fill bark array
          srhoftemp = trhomicro(i)*PI* &
            ((sstemp+barkthick)**2.-sstemp**2.)/(dy*dx)
          tsizescale(ift_index+2,xcor,ycor,k) = &
            (trhof(ift_index+2,xcor,ycor,k)* &
            tsizescale(ift_index+2,xcor,ycor,k)+srhoftemp*barkthick)/ &
            (trhof(ift_index+2,xcor,ycor,k)+srhoftemp)
          tmoist(ift_index+2,xcor,ycor,k) = &
            (trhof(ift_index+2,xcor,ycor,k)* &
            tmoist(ift_index+2,xcor,ycor,k)+srhoftemp*tbarkmoist(i))/ &
            (trhof(ift_index+2,xcor,ycor,k)+srhoftemp)
          trhof(ift_index+2,xcor,ycor,k) = &
            trhof(ift_index+2,xcor,ycor,k)+srhoftemp
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
            do iii=1,10
              do jjj=1,10
                do kkk=1,10
                  test_height= zheight(ii_real,jj_real,kk)+ &
                    (2.*kkk-1.)/20.*(zheight(ii_real,jj_real,kk+1)- &
                    zheight(ii_real,jj_real,kk))
                  xloc = ((ii-1)+(2.*iii-1.)/20.)*dx
                  yloc = ((jj-1)+(2.*jjj-1.)/20.)*dy
                  bot_height = paraboloid(abot,xloc,tlocation(i,j,1), &
                    yloc,tlocation(i,j,2),canopybot)
                  top_height = paraboloid(atop,xloc,tlocation(i,j,1), &
                    yloc,tlocation(i,j,2),canopytop)
                  if (test_height.ge.bot_height.and. &
                    test_height.le.top_height) then 
                    do ift=1,tfuelbins
                      rhoftemp(ift) = rhoftemp(ift)+ &
                        3./2000.*t2bulkdensity(j,ift,i)* &
                        (test_height-canopybot+4.* &
                        (canopytop-canopymaxh)* &
                        ((xloc-tlocation(i,j,1))**2.+ &
                        (yloc-tlocation(i,j,2))**2.)/ &
                        canopydiameter**2.)/(canopytop-canopybot) ! Contribution of one subcell to overall bulk density; 
                    enddo
                  endif
                enddo
              enddo
            enddo
  
            !----- Fill in the 3D arrays for a tree
            do ift=1,tfuelbins
              if (rhoftemp(ift).gt.0) then
                ift_index = (i-1)*ntreefueltypes+ift
                tsizescale(ift_index,ii_real,jj_real,kk) = &
                  (trhof(ift_index,ii_real,jj_real,kk)* &
                  tsizescale(ift_index,ii_real,jj_real,kk)+ &
                  rhoftemp(ift)*t2ss(j,ift,i))/ &
                  (trhof(ift_index,ii_real,jj_real,kk)+rhoftemp(ift))
                tmoist(ift_index,ii_real,jj_real,kk) = &
                  (trhof(ift_index,ii_real,jj_real,kk)* &
                  tmoist(ift_index,ii_real,jj_real,kk)+ &
                  rhoftemp(ift)*t2moisture(j,ift,i))/ &
                  (trhof(ift_index,ii_real,jj_real,kk)+rhoftemp(ift))
                trhof(ift_index,ii_real,jj_real,kk) = &
                  trhof(ift_index,ii_real,jj_real,kk)+rhoftemp(ift)
              endif
            enddo
          enddo
        enddo
      enddo
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
        target_mass = target_mass + &
          t2bulkdensity(j,k,i)*PI*tcrowndiameter(j,i)**2./ &
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
          actual_mass = actual_mass+ &
            trhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
        enddo
      enddo
    enddo
  enddo
  print*,'Trees actual fuel mass:',actual_mass
  print*,'Trees error:',100-actual_mass/target_mass*100,'%'

end subroutine tree_fuels_create 

!-----------------------------------------------------------------------
! litter_fuels_create is a function which computes the characteristics
! of a base litter from trees and fills rhof, sizescale, moisture,
! and fueldepth arrays for litter and subracts mass of the grass 
! arrays for tree shading.
!-----------------------------------------------------------------------
subroutine litter_fuels_create
  use constant_variables, only : pi
  use grid_variables, only : nx,ny,nz,dx,dy,zs,zheight,rhof
  use fuel_variables, only : ntspecies,tfuelbins,grassconstant, &
    litterconstant,trhof,tfuelbins,lrhof,lsizescale,lmoist,lfueldepth, &
    ldepth,lrho,lss,lmoisture,tcrowndiameter,igrass,ngrass,theight, &
    gdepth,ntrees
  use infile_variables, only : infuel
  implicit none
  
  ! Local variables
  integer :: i,j,k,ift,ift_grass
  integer :: zmax
  real :: target_mass,actual_mass
  real :: rhocolumn,coverfactor,shadefactor
  
  ! Executable code
  !----- Place litter on ground and remove grass to account for shading
  print*,'Placing litter and removing grass to account for shading'
  zmax=0
  do ift = 1,ntspecies*tfuelbins
    do i=1,nx
      do j=1,ny
        ! Determine factors for placing litter and removing grass
        rhocolumn = 0
        do k=1,zmax
          rhocolumn = rhocolumn+ &
            sum(trhof((ift-1)*tfuelbins+1:ift*tfuelbins,i,j,k))* &
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
                  rhof(ift_grass+infuel,i,j,k) = &
                    rhof(ift_grass+infuel,i,j,k)*shadefactor
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
              lrhof(ift,i,j,k) = lrho(ift)* &
                (lfueldepth(ift,i,j)+zs(i,j)-zheight(i,j,k))/ &
                (zheight(i,j,k+1)-zheight(i,j,k))
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
    target_mass = target_mass + &
      ntrees(i)*ldepth(i)*PI*tcrowndiameter(1,i)**2.*lrho(i)/4.
  enddo
  print*,'Litter target fuel mass:',target_mass
  actual_mass = 0
  do ift=1,ntspecies
    do i=1,nx
      do j=1,ny
        do k=1,zmax
          actual_mass = actual_mass+ &
            lrhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
        enddo
      enddo
    enddo
  enddo
  print*,'Litter actual fuel mass:',actual_mass
  print*,'Litter error:',100-actual_mass/target_mass*100,'%'

end subroutine litter_fuels_create 
