!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! define_variables definess all variables, both constant and 
! user-defined, used throughout the program
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine define_constant_variables
      !-----------------------------------------------------------------
      ! Constant variables
      !-----------------------------------------------------------------
      use constant_variables
      implicit none

      PI = 3.14159

      end subroutine define_constant_variables

      subroutine define_grid_variables
      !-----------------------------------------------------------------
      ! Grid variables
      !-----------------------------------------------------------------
      use grid_variables
      use infile_variables
      use baseline_variables
      
      implicit none

      integer i,j,k
      integer ii,jj
      integer xbot,xtop,ybot,ytop
      real    cells,xfrac,yfrac
      real,external:: zcart
      real,dimension(2):: xcor,ycor
           
      nfuel = infuel+ngrass+ntspecies*(ilitter+tfuelbins)
      allocate(rhof(nfuel,nx,ny,nz))
      allocate(sizescale(nfuel,nx,ny,nz))
      allocate(moist(nfuel,nx,ny,nz))
      allocate(fueldepth(nfuel,nx,ny,nz))
      
      !-----------------------------------------------------------------
      ! Create topo layer (Should be adjusted for non-flat topo)
      !-----------------------------------------------------------------
      
      allocate(zs(nx,ny))
      allocate(izs(inx,iny))
      allocate(zheight(nx,ny,nz))
      allocate(izheight(inx,iny,inz))
      
      if(ifuelin.eq.1.and.(inx.ne.nx.or.idx.ne.dx.or.
     &  iny.ne.ny.or.idy.ne.dy.or.inz.ne.nz.or.idz.ne.dz
     &  .or.aa1.ne.iaa1)) iintpr=1

      if (topofile.eq.'flat'.or.topofile.eq.'') then ! No topo
        zs(:,:)=0.0
        izs(:,:)=0.0
        print *,'Not using topo'      
      elseif (iintpr.eq.1) then ! Topo with existing fuels
        print *,'Reading given fuel topo file = ',topofile
        open (1,file=topofile,form='unformatted',status='old')
        read (1) izs
        close (1)
        do i=1,nx
          do j=1,ny
            xcor(1) = (i-1)*dx ! Real x lower edge
            xcor(2) = i*dx     ! Real x upper edge
            xbot    = floor(xcor(1)/idx+1) ! Fuel readin grid x lower edge
            xtop    = floor(xcor(2)/idx+1) ! Fuel readin grid x upper edge
            ycor(1) = (j-1)*dy ! Real y lower edge
            ycor(2) = j*dy     ! Real y upper edge
            ybot    = floor(ycor(1)/idy+1) ! Fuel readin grid y lower edge
            ytop    = floor(ycor(2)/idy+1) ! Fuel readin grid y upper edge
            cells   = 0.
            do ii=xbot,xtop
              do jj=ybot,ytop
                xfrac  = (min(ii*idx,xcor(2))-max((ii-1)*idx,xcor(1)))/idx
                yfrac  = (min(jj*idy,ycor(2))-max((jj-1)*idy,ycor(1)))/idy
                cells  = cells+xfrac*yfrac
                zs(i,j)= zs(i,j)+xfrac*yfrac*izs(ii,jj)
              enddo
            enddo
            zs(i,j) = zs(i,j)/cells
          enddo
        enddo
      else ! Normal topo
        print *,'Reading topo file = ',topofile
        open (1,file=topofile,form='unformatted',status='old')
        read (1) zs
        close (1)
      endif
      if(minval(zs).gt.0)then ! Reduce topo values to least common value
        izs = izs-minval(zs)
        zs  = zs-minval(zs)
        open (2,file='toporeduced.dat',form='unformatted',status='unknown')
        write(2) zs
        close(2)
      endif

      do i=1,nx
        do j=1,ny
          do k=1,nz
            if (aa1.eq.0) then
              zheight(i,j,k) = zs(i,j)+(k-1)*dz
            else
              zheight(i,j,k) = zcart(aa1,(k-1)*dz,nz,dz,zs(i,j))
            endif
            if(i.eq.1.and.j.eq.1) print*,'cell',k,'bottom height',zheight(i,j,k)
          enddo
        enddo
      enddo
      if (iintpr.eq.1) then ! Topo with existing fuels
        print*,'Readin fuel grid heights'
        do i=1,inx
          do j=1,iny
            do k=1,inz
              if (iaa1.eq.0) then
                izheight(i,j,k) = izs(i,j)+(k-1)*idz
              else
                izheight(i,j,k) = zcart(iaa1,(k-1)*idz,inz,idz,izs(i,j))
              endif
              if(i.eq.1.and.j.eq.1) print*,'cell',k,'bottom height',izheight(i,j,k)
            enddo
          enddo
        enddo
      endif

      end subroutine define_grid_variables

      subroutine define_baseline_variables
      !-----------------------------------------------------------------
      ! Variables unique to the baseline
      !-----------------------------------------------------------------
      use grid_variables
      use baseline_variables

      implicit none

      integer i,j,ift,itree
      integer,allocatable:: numarray(:)
      real,dimension(7+3*tfuelbins):: temp_array
      
      if (igrass.ne.0) then
        allocate(grhof(ngrass,nx,ny,nz))
        allocate(gsizescale(ngrass,nx,ny,nz))
        allocate(gmoist(ngrass,nx,ny,nz))
        allocate(gfueldepth(ngrass,nx,ny))
      endif
      if (itrees.ne.0) then
        allocate(trhof(ntspecies*tfuelbins,nx,ny,nz))
        allocate(tsizescale(ntspecies*tfuelbins,nx,ny,nz))
        allocate(tmoist(ntspecies*tfuelbins,nx,ny,nz))
        allocate(tfueldepth(ntspecies*tfuelbins,nx,ny))
      endif
      if (ilitter.ne.0) then
        allocate(lrhof(ntspecies,nx,ny,nz))
        allocate(lsizescale(ntspecies,nx,ny,nz))
        allocate(lmoist(ntspecies,nx,ny,nz))
        allocate(lfueldepth(ntspecies,nx,ny))
      endif

      !-----------------------------------------------------------------
      ! Groundfuel variables unique to the ground fuels baseline
      !-----------------------------------------------------------------
      if (igrass.eq.1) then
        allocate(grho(ngrass))
        allocate(gmoisture(ngrass))
        allocate(gss(ngrass))
        allocate(gdepth(ngrass))
        
        open (1,file=grassfile)
        read (1,*) grho       ! bulk density of grass [kg/m3]
        read (1,*) gmoisture  ! moisture content of grass
        read (1,*) gss        ! size scale of grass [m]
        read (1,*) gdepth     ! depth of grass [m]
        close (1)
      elseif (igrass.eq.2) then
        print*,'Reading 2D ground fuel file'
        open (1,file=grassfile,form='unformatted',status='old')
        do ift=1,ngrass
          read (1) grhof(ift,:,:,1)      ! bulk density of grass [kg/m3]
          read (1) gmoist(ift,:,:,1)     ! moisture content of grass
          read (1) gsizescale(ift,:,:,1) ! size scale of grass [m]
          read (1) gfueldepth(ift,:,:) ! depth of grass [m]
        enddo
        close (1)
      endif
      
      !-----------------------------------------------------------------
      ! Tree variables unique to the tree baseline
      !-----------------------------------------------------------------
      if (itrees.eq.1) then
        allocate(ntrees(ntspecies)) ! Number of trees for each species
        allocate(tstemdensity(ntspecies)) ! Stem density of each species [stem/ha]
        allocate(theight(2,ntspecies)) ! Tree heights [m]
        allocate(tcrownbotheight(2,ntspecies)) ! Height to live crown [m]
        allocate(tcrowndiameter(2,ntspecies)) ! Crown diameter [m]
        allocate(tcrownmaxheight(2,ntspecies)) ! Height to max crown diameter [m]
        allocate(t1bulkdensity(tfuelbins,ntspecies)) ! Crown fuel bulk density [kg/m3]
        allocate(t1moisture(tfuelbins,ntspecies)) ! Crown fuel moisture content [fraction]
        allocate(t1ss(tfuelbins,ntspecies)) ! Crown fuel size scale [m]
      
        open (2,file=treefile)
        read (2,*) tstemdensity
        read (2,*) theight(1,:)
        read (2,*) theight(2,:)
        read (2,*) tcrownbotheight(1,:)
        read (2,*) tcrownbotheight(2,:)
        read (2,*) tcrowndiameter(1,:)
        read (2,*) tcrowndiameter(2,:)
        read (2,*) tcrownmaxheight(1,:)
        read (2,*) tcrownmaxheight(2,:)
        do i=1,tfuelbins
          read(2,*) t1bulkdensity(i,:)
          read(2,*) t1moisture(i,:)
          read(2,*) t1ss(i,:)
        enddo
        close (2)
      endif
      
      if (itrees.eq.2.or.itrees.eq.3) then
        itree = 0
        open (2,file=treefile)
        do
          read (2,*,end=10)
          itree = itree+1
        enddo
10      rewind(2)
        
        allocate(tspecies(itree))
        do i=1,itree
          read(2,*) temp_array(:)
          tspecies(i)=temp_array(1)
        enddo
        rewind(2)
        ntspecies = maxval(tspecies)
        allocate(ntrees(ntspecies)) ! Total number of trees for each species
        ntrees=0
        do i=1,itree
          ntrees(tspecies(i)) = ntrees(tspecies(i))+1
        enddo
        allocate(tlocation(ntspecies,maxval(ntrees),2)) ! Tree cartesian coordinates [m,m]
        allocate(theight(maxval(ntrees),ntspecies)) ! Tree heights [m]
        allocate(tcrownbotheight(maxval(ntrees),ntspecies)) ! Height to live crown [m]
        allocate(tcrowndiameter(maxval(ntrees),ntspecies)) ! Crown diameter [m]
        allocate(tcrownmaxheight(maxval(ntrees),ntspecies)) ! Height to max crown diameter [m]
        allocate(t2bulkdensity(maxval(ntrees),tfuelbins,ntspecies)) ! Crown fuel bulk density [kg/m3]
        allocate(t2moisture(maxval(ntrees),tfuelbins,ntspecies)) ! Crown fuel moisture content [fraction]
        allocate(t2ss(maxval(ntrees),tfuelbins,ntspecies)) ! Crown fuel size scale [m]
        allocate(numarray(ntspecies))
        numarray(:)=0
        do i=1,itree
          read(2,*) temp_array(:)
          numarray(tspecies(i)) = numarray(tspecies(i))+1
          tlocation(tspecies(i),numarray(tspecies(i)),:) = temp_array(2:3)
          theight(numarray(tspecies(i)),tspecies(i)) = temp_array(4)
          tcrownbotheight(numarray(tspecies(i)),tspecies(i)) = temp_array(5)
          tcrowndiameter(numarray(tspecies(i)),tspecies(i)) = temp_array(6)
          tcrownmaxheight(numarray(tspecies(i)),tspecies(i)) = temp_array(7)
          do j=1,tfuelbins
            t2bulkdensity(numarray(tspecies(i)),j,tspecies(i)) = temp_array(8+3*(j-1))
            t2moisture(numarray(tspecies(i)),j,tspecies(i)) = temp_array(9+3*(j-1))
            t2ss(numarray(tspecies(i)),j,tspecies(i)) = temp_array(10+3*(j-1))
          enddo
        enddo
        close (2)
      endif

      if (ilitter.eq.1) then
        allocate(lrho(ntspecies))
        allocate(lmoisture(ntspecies))
        allocate(lss(ntspecies))
        allocate(ldepth(ntspecies))
      
        open (3,file=litterfile)
        read (3,*) lrho       ! bulk density of litter [kg/m3]
        read (3,*) lmoisture  ! moisture content of litter [fraction]
        read (3,*) lss        ! size scale of litter [m]
        read (3,*) ldepth     ! depth of litter [m]
        close (3)
      endif

      end subroutine define_baseline_variables
