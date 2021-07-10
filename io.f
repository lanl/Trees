!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! io contains the functions for reading input namelists and tree data 
! files and writing .dat files
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine namelist_input
      !-----------------------------------------------------------------
      ! namelist_input is a function reads in the user-defined variables 
      ! from the treelist then assigns those variables for use by the 
      ! trees softeware.
      !-----------------------------------------------------------------
      use grid_variables
      use infile_variables
      use baseline_variables
      use treatment_variables
      implicit none

      namelist/fuellist/
     .   nx,ny,nz,dx,dy,dz,aa1,topofile,
     .   ifuelin,rhoffile,moistfile,ssfile,afdfile,
     .   inx,iny,inz,idx,idy,idz,iaa1,infuel,
     .   igrass,ngrass,grassconstant,grassfile,
     .   itrees,ntspecies,tfuelbins,tdnx,tdny,treefile,
     .   ilitter,litterconstant,litterfile,
     .   itreatment,sdnx,sdny,
     .   sdiameter,sheight,sprho,smoist,sdepth
      
      ! Area of interest arrays need to be allocated before calling namelist
      allocate(tdnx(2)) ! Array of the cell range (x)  where the trees are applied
      allocate(tdny(2)) ! Array of the cell range (x)  where the trees are applied
      allocate(sdnx(2)) ! Array of the cell range (x)  where the treatment is applied
      allocate(sdny(2)) ! Array of the cell range (x)  where the treatment is applied
      
      ! Set the default values that will be overwritten by the namelist if present
      aa1 = 0.1
      iaa1= -1
      topofile = 'flat'
      grassconstant = 5
      ntspecies = 1
      tfuelbins = 1
      litterconstant = 5

      open(unit=15,file='fuellist',form='formatted',status='old')
           read (15,nml=fuellist)
      close(15)

      ! Corrections for if variables not specifiedi on namelist
      if (tdnx(1).eq.0) then
        tdnx(1) = 1
        tdnx(2) = dx*nx
        tdny(1) = 1
        tdny(2) = dy*ny
      endif
      if (sdnx(1).eq.0) then
        sdnx(1) = 1
        sdnx(2) = dx*nx
        sdny(1) = 1
        sdny(2) = dy*ny
      endif
      tdnx(1)=floor(tdnx(1)/dx+1)
      tdnx(2)=ceiling(tdnx(2)/dx)
      tdny(1)=floor(tdny(1)/dy+1)
      tdny(2)=ceiling(tdny(2)/dy)
      sdnx(1)=floor(sdnx(1)/dx+1)
      sdnx(2)=ceiling(sdnx(2)/dx)
      sdny(1)=floor(sdny(1)/dy+1)
      sdny(2)=ceiling(sdny(2)/dy)

      if (ifuelin.eq.1) then
        if(inx.eq.0) inx=nx
        if(iny.eq.0) iny=ny
        if(inz.eq.0) inz=nz
        if(idx.eq.0) idx=dx
        if(idy.eq.0) idy=dy
        if(idz.eq.0) idz=dz
        if(iaa1.eq.-1) iaa1=aa1
        if(infuel.eq.0) infuel=1
      endif

      end subroutine namelist_input

      subroutine output
      !-----------------------------------------------------------------
      ! output is a function which writes the .dat files for use in 
      ! FIRETEC or QUIC-Fire
      !-----------------------------------------------------------------
      use grid_variables
      implicit none
      integer ift,i,j,k,lfuel
      real,dimension(nfuel):: nonzero
      
      nonzero(:) = 1
      lfuel = 1
      do ift=1,nfuel
        if (sum(rhof(ift,:,:,:)).le.0) then
          nonzero(ift) = 0
        else
          do k=1,nz
            if (sum(rhof(ift,:,:,k)).gt.0) lfuel = max(lfuel,k)
          enddo
        endif
      enddo

      print*,'Exporting data to .dat files'
      open (1,file='treesrhof.dat',form='unformatted',status='unknown')
      do ift=1,nfuel
        if (nonzero(ift).ne.0)  write (1) rhof(ift,:,:,:)
      enddo
      close (1)

      open (1,file='treesfueldepth.dat',form='unformatted',status='unknown')
      do ift=1,nfuel
        if (nonzero(ift).ne.0)  write (1) fueldepth(ift,:,:,:)
      enddo
      close (1)

      open (1,file='treesss.dat',form='unformatted',status='unknown')
      do ift=1,nfuel
        if (nonzero(ift).ne.0)  write (1) sizescale(ift,:,:,:)
      enddo
      close (1)

      open (1,file='treesmoist.dat',form='unformatted',status='unknown')
      do ift=1,nfuel
        if (nonzero(ift).ne.0)  write (1) moist(ift,:,:,:)
      enddo
      close (1)

      print*,'Your nfuel is',int(sum(nonzero(:)))
      print*,'Your lfuel is',lfuel
      
      end subroutine output
      
      subroutine fuel_readin
      !-----------------------------------------------------------------
      ! fuel_readin is a function which reads in .dat files for use in 
      ! FIRETEC or QUIC-Fire
      !-----------------------------------------------------------------
      use grid_variables
      use infile_variables
      use baseline_variables
      implicit none
      integer ift,i,j,k
      integer ii,jj,kk,kkk
      integer xbot,xtop,ybot,ytop,zbot,ztop
      real cells,xfrac,yfrac,zfrac
      real zstemp,rhoftemp 
      real target_mass,actual_mass
      real,dimension(2):: xcor,ycor,zcor
      
      allocate(irhof(infuel,inx,iny,inz))
      allocate(imoist(infuel,inx,iny,inz))
      allocate(iss(infuel,inx,iny,inz))
      allocate(iafd(infuel,inx,iny,inz))

      ! Read in fuel files
      print*,"Reading existing fuel files"
      open(unit=1,file=rhoffile,form='unformatted',status='unknown')
        read(1) irhof
      close(1)
      where (irhof<0)
        irhof = 0
      endwhere
      print*,'irhof','min=',minval(irhof),'max=',maxval(irhof),'avg=',sum(irhof)/(inx*iny*inz)
      open(unit=2,file=moistfile,form='unformatted',status='unknown')
        read(2) imoist
      close(2)
      where (imoist<0)
        imoist = 0
      endwhere
      print*,'imoist','min=',minval(imoist),'max=',maxval(imoist),'avg=',sum(imoist)/(inx*iny*inz)
      open(unit=3,file=ssfile,form='unformatted',status='unknown')
        read(3) iss
      close(3)
      if (ssfile.eq."sav.dat.orig")then ! Special case if ss file was sav
        do i=1,inx
          do j=1,iny
            do k=1,inz
              do ift=1,infuel
                if(iss(ift,i,j,k).gt.1) iss(ift,i,j,k)=2./iss(ift,i,j,k)
              enddo
            enddo
          enddo
        enddo
      endif
      where (iss<0)
        iss = 0
      endwhere
      print*,'iss','min=',minval(iss),'max=',maxval(iss),'avg=',sum(iss)/(inx*iny*inz)
      open(unit=4,file=afdfile,form='unformatted',status='unknown')
        read(4) iafd
      close(4)
      where (iafd<0)
        iafd = 0
      endwhere
      print*,'iafd','min=',minval(iafd),'max=',maxval(iafd),'avg=',sum(iafd)/(inx*iny*inz)

      print*,"Exisiting fuel files readin"

      ! Interpolate read file onto FIRETEC grid
      if (iintpr.eq.1) then
        print*,"Interpolating readin fuel files to desired FIRETEC grid"
        do i=1,nx
          do j=1,ny
            xcor(1) = (i-1)*dx ! Real x lower edge
            xcor(2) = i*dx     ! Real x upper edge
            xbot    = floor(xcor(1)/idx+1) ! Fuel readin grid x lower edge
            xtop    = min(inx,floor(xcor(2)/idx+1)) ! Fuel readin grid x upper edge
            ycor(1) = (j-1)*dy ! Real y lower edge
            ycor(2) = j*dy     ! Real y upper edge
            ybot    = floor(ycor(1)/idy+1) ! Fuel readin grid y lower edge
            ytop    = min(iny,floor(ycor(2)/idy+1)) ! Fuel readin grid y upper edge
            do k=1,nz-1
              cells = 0.
              do ii=xbot,xtop
                do jj=ybot,ytop
                  do kk=1,inz-1
                    if(izheight(ii,jj,kk).ge.zheight(i,j,k+1)) exit
                    if(izheight(ii,jj,kk).ge.zheight(i,j,k))then
                      do ift=1,infuel
                        xfrac = (min(ii*idx,xcor(2))-max((ii-1)*idx,xcor(1)))/idx
                        yfrac = (min(jj*idy,ycor(2))-max((jj-1)*idy,ycor(1)))/idy
                        zfrac = (min(izheight(ii,jj,kk+1),zheight(i,j,k+1))-max(izheight(ii,jj,kk),zheight(i,j,k)))/
     +                    (izheight(ii,jj,kk+1)-izheight(ii,jj,kk))
                        cells = cells+xfrac*yfrac*zfrac
                        rhoftemp=irhof(ift,ii,jj,kk)*xfrac*yfrac*zfrac
                        if(rhoftemp.gt.0)then
                        sizescale(ift,i,j,k) = (rhof(ift,i,j,k)*sizescale(ift,i,j,k)+rhoftemp*iss(ift,ii,jj,kk))/
     +                    (rhof(ift,i,j,k)+rhoftemp)
                        moist(ift,i,j,k) = (rhof(ift,i,j,k)*moist(ift,i,j,k)+rhoftemp*imoist(ift,ii,jj,kk))/
     +                    (rhof(ift,i,j,k)+rhoftemp)
                        fueldepth(ift,i,j,k) = (rhof(ift,i,j,k)*fueldepth(ift,i,j,k)+rhoftemp*iafd(ift,ii,jj,kk))/
     +                    (rhof(ift,i,j,k)+rhoftemp)
                        rhof(ift,i,j,k)=rhof(ift,i,j,k)+rhoftemp
                        endif
                      enddo
                    endif
                  enddo
                enddo
              enddo
              do ift=1,infuel
                if(cells.gt.0) rhof(ift,i,j,k)=rhof(ift,i,j,k)/cells
              enddo
            enddo
          enddo
        enddo
      endif
      
      ! Print out the target and actual fuel masses for comparisons sake
      target_mass = 0
      do ift=1,infuel
        do i=1,min(inx,nint(nx/dx*idx))
          do j=1,min(iny,nint(ny/dy*idy))
            do k=1,inz-1
              target_mass = target_mass+irhof(ift,i,j,k)*idx*idy*(izheight(i,j,k+1)-izheight(i,j,k))
            enddo
          enddo
        enddo
      enddo
      print*,'Readin target fuel mass:',target_mass
      actual_mass = 0
      do ift=1,infuel
        do i=1,nx
          do j=1,ny
            do k=1,nz-1
              actual_mass = actual_mass+rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
            enddo
          enddo
        enddo
      enddo
      print*,'Readin actual fuel mass:',actual_mass
      print*,'Readin error:',actual_mass/target_mass*100,'%'
       
      deallocate(irhof)
      deallocate(imoist)
      deallocate(iss)
      deallocate(iafd)
      
      end subroutine fuel_readin
