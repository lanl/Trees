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
     .   nx,ny,nz,dx,dy,dz,aa1,singlefuel,topofile,
     .   ifuelin,rhoffile,moistfile,ssfile,afdfile,
     .   inx,iny,inz,idx,idy,idz,iaa1,infuel,
     .   igrass,ngrass,grassconstant,grassfile,
     .   itrees,ntspecies,tfuelbins,tdnx,tdny,treefile,istem,
     .   ndatax,ndatay,datalocx,datalocy, !JSM added for populate function
     .   ilitter,litterconstant,litterfile,
     .   itreatment,sdnx,sdny,
     .   sdiameter,sheight,sprho,smoist,sdepth
      
      ! Area of interest arrays need to be allocated before calling namelist
      allocate(tdnx(2)); tdnx(:)=0 ! Array of the cell range (x)  where the trees are applied
      allocate(tdny(2)); tdny(:)=0 ! Array of the cell range (x)  where the trees are applied
      allocate(sdnx(2)); sdnx(:)=0 ! Array of the cell range (x)  where the treatment is applied
      allocate(sdny(2)); sdny(:)=0 ! Array of the cell range (x)  where the treatment is applied
      
      ! Set the default values that will be overwritten by the namelist if present
      aa1           =0.1
      singlefuel    =0.0
      topofile      ='flat'
      ifuelin       =0
      iaa1          =-1
      igrass        =0
      ngrass        =1
      grassconstant =5
      itrees        =0
      ntspecies     =1
      tfuelbins     =1
      istem         =0
      ndatax        =0
      ndatay        =0
      datalocx      =0  !JSM added for populate function
      datalocy      =0  !JSM added for populate function
      ilitter       =0
      litterconstant=5
      itreatment    =0

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
      
      if (ndatax.eq.0) ndatax=nx*dx
      if (ndatay.eq.0) ndatay=ny*dy

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

      !Find number species if using FastFuels dataset
      if (itrees.eq.7) then
        call find_fastfuels_numspecies
        tfuelbins = 1
      endif

      end subroutine namelist_input

      subroutine output_nfuel
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
      
      end subroutine output_nfuel
      
      subroutine output_1fuel
      !-----------------------------------------------------------------
      ! output_1fuel is a function which writes the .dat files for use in 
      ! FIRETEC or QUIC-Fire but combines all fuels into a single fuel
      ! type.
      !-----------------------------------------------------------------
      use grid_variables
      implicit none
      integer ift,i,j,k,lfuel
      real,allocatable :: srhof(:,:,:),sss(:,:,:),smoist(:,:,:),safd(:,:,:)
      
      allocate(srhof(nx,ny,nz))
      allocate(sss(nx,ny,nz))
      allocate(smoist(nx,ny,nz))
      allocate(safd(nx,ny,nz))
      do i=1,nx
        do j=1,ny
          do k=1,nz
            do ift=1,nfuel
              if(rhof(ift,i,j,k).gt.0)then
                sss(i,j,k)=(sss(i,j,k)*srhof(i,j,k)+sizescale(ift,i,j,k)*rhof(ift,i,j,k))
     &            /(srhof(i,j,k)+rhof(ift,i,j,k))
                smoist(i,j,k)=(smoist(i,j,k)*srhof(i,j,k)+moist(ift,i,j,k)*rhof(ift,i,j,k))
     &            /(srhof(i,j,k)+rhof(ift,i,j,k))
                safd(i,j,k)=(safd(i,j,k)*srhof(i,j,k)+fueldepth(ift,i,j,k)*rhof(ift,i,j,k))
     &            /(srhof(i,j,k)+rhof(ift,i,j,k))
                srhof(i,j,k)=srhof(i,j,k)+rhof(ift,i,j,k)
              endif
            enddo
          enddo
        enddo
      enddo

      print*,'Exporting data to .dat files'
      open (1,file='treesrhof.dat',form='unformatted',status='unknown')
      write (1) srhof
      close (1)

      open (1,file='treesfueldepth.dat',form='unformatted',status='unknown')
      write (1) safd
      close (1)

      open (1,file='treesss.dat',form='unformatted',status='unknown')
      write (1) sss
      close (1)

      open (1,file='treesmoist.dat',form='unformatted',status='unknown')
      write (1) smoist
      close (1)

      print*,'Your nfuel is 1'
      
      end subroutine output_1fuel


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine will find the number of species in the fast fuels data
! Added 10/18 JO
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine find_fastfuels_numspecies
        !-----------------------------------------------------------------
        !-----------------------------------------------------------------
            use grid_variables
            use baseline_variables      
            implicit none

            !variables for finding # of fuel speies
            integer i,j,ift,ff_len,min_val_sp, max_val_sp
            integer,allocatable :: final_uni_sp(:), uni_sp(:)
            real,dimension(19):: temp_array ! FF trees csv has at least 19 columns

            ff_len = 0
            open (2,file=treefile)
            do
               read (2,*,end=19) !length of FF columns
               ff_len = ff_len+1
            enddo
19          rewind(2)

            allocate(tspecies(ff_len-1))
            read(2,*) !read 1st line and throw away, has column headers
            do i=1,ff_len-1
               read(2,*) temp_array(:)
               tspecies(i)=temp_array(5) !take from sp_grp, 5th pos
               !print*,'species# = ',temp_array(5)
            enddo
            rewind(2)
      
            !find unique number of tree species!
            allocate(uni_sp(maxval(tspecies)))
            min_val_sp = minval(tspecies)-1
            max_val_sp = maxval(tspecies)
            i=0
            do while (min_val_sp<max_val_sp)
                i = i+1
                min_val_sp = minval(tspecies, mask=tspecies>min_val_sp)
                uni_sp(i) = min_val_sp
            enddo
            allocate(final_uni_sp(i), source=uni_sp(1:i)) 
            ntspecies = count(final_uni_sp==final_uni_sp)
            print*,'New Set Species = ',ntspecies

            deallocate(tspecies)
            deallocate(uni_sp)
            deallocate(final_uni_sp)


        end subroutine find_fastfuels_numspecies
      
