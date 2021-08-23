!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! defines subroutines for reading and interpreting various types of 
! tree data files
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine treesGeneral_readin
      !-----------------------------------------------------------------
      ! treesGeneral_readin is a function which reads in a general trees
      ! file for use in FIRETEC or QUIC-Fire
      !-----------------------------------------------------------------
      use grid_variables
      use baseline_variables
      
      implicit none
      integer i
        
      allocate(ntrees(ntspecies)); ntrees(:)=0.0 ! Number of trees for each species
      allocate(tstemdensity(ntspecies)) ! Stem density of each species [stem/ha]
      allocate(theight(2,ntspecies)) ! Tree heights [m]
      allocate(tcrownbotheight(2,ntspecies)) ! Height to live crown [m]
      allocate(tcrowndiameter(2,ntspecies)) ! Crown diameter [m]
      allocate(tcrownmaxheight(2,ntspecies)) ! Height to max crown diameter [m]
      allocate(t1bulkdensity(tfuelbins,ntspecies)) ! Crown fuel bulk density [kg/m3]
      allocate(t1moisture(tfuelbins,ntspecies)) ! Crown fuel moisture content [fraction]
      allocate(t1ss(tfuelbins,ntspecies)) ! Crown fuel size scale [m]
      if(istem.eq.1)then
        allocate(trhomicro(ntspecies)) ! Micro-density of tree species [kg/m3]
        allocate(tdbh(2,ntspecies)) ! Tree diameter breast heights [m]
        allocate(tstemmoist(ntspecies)) ! Tree stem moisture content [fraction]
        allocate(tbarkthick(2,ntspecies)) ! Bark Thickness [m]
        allocate(tbarkmoist(ntspecies)) ! Bark moisture content [fraction]
      endif

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
      if(istem.eq.1)then
        read (2,*) trhomicro(:)
        read (2,*) tdbh(1,:)
        read (2,*) tdbh(2,:)
        read (2,*) tstemmoist(:)
        read (2,*) tbarkthick(1,:)
        read (2,*) tbarkthick(2,:)
        read (2,*) tbarkmoist(:)
      endif 
      do i=1,tfuelbins
        read(2,*) t1bulkdensity(i,:)
        read(2,*) t1moisture(i,:)
        read(2,*) t1ss(i,:)
      enddo
      close (2)
      
      end subroutine treesGeneral_readin
      
      subroutine treelist_readin
      !-----------------------------------------------------------------
      ! treesGeneral_readin is a function which reads in a treelist
      ! file for use in FIRETEC or QUIC-Fire
      !-----------------------------------------------------------------
      use grid_variables
      use baseline_variables
      
      implicit none
      integer i,j,ift,itree
      integer,allocatable:: numarray(:)
      real,dimension(7+3*tfuelbins):: temp_array
      
      !!!------JSM: Parameters/Variables for populate function-----!!!
      real:: nsub,nsubdecimal,rnumx,rnumy,newx,newy
      integer:: q,t,r,s,tindex,dataleft,dataright,databottom,datatop,treecount,num
      integer,allocatable:: rounddown(:),ntreesold(:)
        
      itree = 0
      open (2,file=treefile)
      do
        read (2,*,end=10)
        itree = itree+1
      enddo
10    rewind(2)
      
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

  !!!!!------JSM added the following block of code for populate------!!!!!

      !---Determine how many dataset subdomains fit within your main domain
      !if(ndatax.lt.nx*dx.or.ndatay.lt.ny*dy) then
      nsub = (nx*dx*ny*dy)/(ndatax*ndatay)
      print*,'Number of subdomains = ',nsub
      !endif
      allocate(ntreesold(ntspecies))
      ntreesold = ntrees

      !---If there is an integer number of subdomains, just multiply the number
      !of trees in each column of ntrees by that number of subdomains (minus 1
      !for the original dataset trees that are already in ntrees) - if it is not
      !an integer, we need to determine how many of the trees of each species
      !that we need to multiply by the integer above the decimal, and those we
      !need to multiply by the integer below the decimal...
      if(nsub.eq.nint(nsub)) then
         ntrees = ntrees*(nsub)
         print*,'ntrees = ',ntrees

      else
         allocate(rounddown(ntspecies)) 
         print*,'Not an integer number of subdomains...'
         nsubdecimal = (nint(nsub)-nsub)
         if(nsubdecimal.lt.0) then
             rounddown = nint((1-abs(nsubdecimal))*ntrees)
         else
             rounddown = nint(nsubdecimal*ntrees)
         endif
         print*,'Number of trees to replicate rounded down for each species',rounddown

         print*,'old ntrees = ',ntrees
         print*,'rounddown integer = ',floor(nsub)
         print*,'roundup integer = ',ceiling(nsub)
         do i=1,ntspecies
             ntrees(i) = rounddown(i)*(floor(nsub)) + (ntreesold(i)-rounddown(i))*(ceiling(nsub))
         enddo
         print*,'new ntrees = ',ntrees
      endif
  !!!!!------END of JSM additions for populate-------!!!!!

      allocate(tlocation(ntspecies,maxval(ntrees),2)); tlocation(:,:,:)=0.0 ! Tree cartesian coordinates [m,m]
      allocate(theight(maxval(ntrees),ntspecies)); theight(:,:)=0.0 ! Tree heights [m]
      allocate(tcrownbotheight(maxval(ntrees),ntspecies)); tcrownbotheight(:,:)=0.0 ! Height to live crown [m]
      allocate(tcrowndiameter(maxval(ntrees),ntspecies)); tcrowndiameter(:,:)=0.0 ! Crown diameter [m]
      allocate(tcrownmaxheight(maxval(ntrees),ntspecies)); tcrownmaxheight(:,:)=0.0 ! Height to max crown diameter [m]
      allocate(t2bulkdensity(maxval(ntrees),tfuelbins,ntspecies)); t2bulkdensity(:,:,:)=0.0 ! Crown fuel bulk density [kg/m3]
      allocate(t2moisture(maxval(ntrees),tfuelbins,ntspecies)); t2moisture(:,:,:)=0.0 ! Crown fuel moisture content [fraction]
      allocate(t2ss(maxval(ntrees),tfuelbins,ntspecies)); t2ss(:,:,:)=0.0 ! Crown fuel size scale [m]
      allocate(numarray(ntspecies)); numarray(:)=0
      do i=1,itree
        read(2,*) temp_array(:)
        numarray(tspecies(i)) = numarray(tspecies(i))+1
        tlocation(tspecies(i),numarray(tspecies(i)),1) = temp_array(2)+datalocx
        tlocation(tspecies(i),numarray(tspecies(i)),2) = temp_array(3)+datalocy
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

      !!!-------JSM More populate stuff-------!!!
      ! For each of the trees in the dataset which are placed above, we copy the
      ! information for each of them into the arrays and choose a new location;
      ! if the initial random location is within the chosen dataset area, we
      ! choose another one

      !First we need to find the area for where the dataset will live

      dataleft = datalocx
      dataright = (datalocx + ndatax)
      databottom = datalocy
      datatop = (datalocy + ndatay)

      print*,'dataset lives at these coordinates: ',dataleft,dataright,databottom,datatop

      if (nsub.gt.1) then
          if (nsub.eq.int(nsub)) then
             do q=1,ntspecies
                tindex = ntreesold(q)
                do r=1,ntreesold(q)
                    do s=1,int(nsub)-1         
                       tindex = tindex+1
                       newx = tlocation(q,r,1)
                       newy = tlocation(q,r,2)
                       !print*,'old location = ',newx,newy   
                       do while (newx.ge.dataleft.and.newx.le.dataright.and.newy.ge.databottom.and.newy.le.datatop.or.newx.gt.nx*dx.or.newx.lt.0.or.newy.gt.ny*dy.or.newy.lt.0)
                          call random_number(rnumx)
                          newx = rnumx*nx*dx
                          call random_number(rnumy)
                          newy = rnumy*ny*dy
                       enddo

                    tlocation(q,tindex,1) = newx
                    tlocation(q,tindex,2) = newy
                    theight(tindex,q) = theight(r,q)
                    tcrownbotheight(tindex,q) = tcrownbotheight(r,q)
                    tcrowndiameter(tindex,q) = tcrowndiameter(r,q)
                    tcrownmaxheight(tindex,q) = tcrownmaxheight(r,q)
                    do j=1,tfuelbins
                       t2bulkdensity(tindex,j,q) = t2bulkdensity(r,j,q)
                       t2moisture(tindex,j,q) = t2moisture(r,j,q)
                       t2ss(tindex,j,q) = t2ss(r,j,q)
                    enddo
                    !print*,'theight old = ',theight(r,q)
                    !print*,'theight new = ',theight(r,q)
                    enddo
                enddo
             enddo
          else
           do q=1,ntspecies
           tindex = ntreesold(q)
             do r=1,rounddown(q)
                do s=1,floor(nsub)-1
                  tindex=tindex+1
                  !print*,'original index = ',r
                  !print*,'tindex = ',tindex
                  !Choose a new location
                  newx = tlocation(q,r,1)
                  newy = tlocation(q,r,2)
                  !print*,'old location = ',newx,newy   
                  do while (newx.ge.dataleft.and.newx.le.dataright.and.newy.ge.databottom.and.newy.le.datatop.or.newx.gt.nx*dx.or.newx.lt.0.or.newy.gt.ny*dy.or.newy.lt.0)
                     call random_number(rnumx)
                     newx = rnumx*nx*dx
                     call random_number(rnumy)
                     newy = rnumy*ny*dy
                  enddo
                 
                  !do while(newy.ge.databottom.and.newy.le.datatop.and.newx.le.dataright.and.newx.ge.dataleft.or.newy.gt.ny*dy.or.newy.lt.0)
                  !   call random_number(rnum)
                  !   newy = rnum*ny*dy
                  !enddo
                  !print*,'new location = ',newx,newy
                  tlocation(q,tindex,1) = newx
                  tlocation(q,tindex,2) = newy
                  theight(tindex,q) = theight(r,q)
                  tcrownbotheight(tindex,q) = tcrownbotheight(r,q)
                  tcrowndiameter(tindex,q) = tcrowndiameter(r,q)
                  tcrownmaxheight(tindex,q) = tcrownmaxheight(r,q)
                  do j=1,tfuelbins
                    t2bulkdensity(tindex,j,q) = t2bulkdensity(r,j,q)
                    t2moisture(tindex,j,q) = t2moisture(r,j,q)
                    t2ss(tindex,j,q) = t2ss(r,j,q)
                  enddo
                  !print*,'theight old = ',theight(r,q)
                  !print*,'theight new = ',theight(r,q)
              enddo
           enddo
          do r=rounddown(q)+1,ntreesold(q)
              do s=1,ceiling(nsub)-1
                  tindex=tindex+1
                  !print*,'original index = ',r
                  !print*,'tindex = ',tindex
                  !Choose a new location
                  newx = tlocation(q,r,1)
                  newy = tlocation(q,r,2)
                  !print*,'old location = ',newx,newy

                  do while (newx.ge.dataleft.and.newx.le.dataright.and.newy.ge.databottom.and.newy.le.datatop.or.newx.gt.nx*dx.or.newx.lt.0.or.newy.gt.ny*dy.or.newy.lt.0)
                     call random_number(rnumx)
                     newx = rnumx*nx*dx
                     call random_number(rnumy)
                     newy = rnumy*ny*dy
                  enddo


                  !print*,'new location = ',newx,newy
                  tlocation(q,tindex,1) = newx
                  tlocation(q,tindex,2) = newy
                  theight(tindex,q) = theight(r,q)
                  tcrownbotheight(tindex,q) = tcrownbotheight(r,q)
                  tcrowndiameter(tindex,q) = tcrowndiameter(r,q)
                  tcrownmaxheight(tindex,q) = tcrownmaxheight(r,q)
                  do j=1,tfuelbins
                    t2bulkdensity(tindex,j,q) = t2bulkdensity(r,j,q)
                    t2moisture(tindex,j,q) = t2moisture(r,j,q)
                    t2ss(tindex,j,q) = t2ss(r,j,q)
                  enddo
                  !print*,'theight old = ',theight(r,q)
                  !print*,'theight new = ',theight(r,q)

              enddo
          enddo
      enddo
      endif
      treecount = 0
      do q=1,ntspecies
          treecount = treecount + ntrees(q)
      enddo
      print*,'Treecount = ',treecount

      num=0
      do q=1,ntspecies
          do t=ntreesold(q)+1,ntrees(q)
              do r=1,ntspecies
                  do s=1,ntrees(r)
                      if (q.ne.r.or.t.ne.s) then
                         do while(abs(tlocation(q,t,1)-tlocation(r,s,1)).lt.0.1.and.abs(tlocation(q,t,2)-tlocation(r,s,2)).lt.0.1)
                            call random_number(rnumx)
                            newx = rnumx*nx*dx
                            call random_number(rnumy)
                            newy = rnumy*ny*dy

                            do while (newx.ge.dataleft.and.newx.le.dataright.and.newy.ge.databottom.and.newy.le.datatop.or.newx.gt.nx*dx.or.newx.lt.0.or.newy.gt.ny*dy.or.newy.lt.0)
                               call random_number(rnumx)
                               newx = rnumx*nx*dx
                               call random_number(rnumy)
                               newy = rnumy*ny*dy
                            enddo

                            num = num+1
                            tlocation(q,t,1) = newx
                            tlocation(q,t,2) = newy
                         enddo
                      endif
                  enddo
               enddo
           enddo
      enddo

      print*,'Number of relocation due to crowding = ',num
      endif
      
      !!!---------END OF JSM ADDITIONS FOR POPULATE----------!!!

      end subroutine treelist_readin
