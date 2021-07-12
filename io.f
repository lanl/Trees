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
