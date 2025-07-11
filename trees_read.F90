!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! defines subroutines for reading and interpreting various types of 
! tree data files
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
! treesGeneral_readin is a function which reads in a treelist
! file for use in FIRETEC or QUIC-Fire
!-----------------------------------------------------------------------
subroutine treelist_readin
  use grid_variables, only : nx,ny,dx,dy,ndatax,ndatay,datalocx,datalocy
  use fuel_variables, only : tfuelbins,newtreefile,ntspecies,ntrees, &
    tlocation,theight,tcrownbotheight,tcrownmaxheight,itrees,treefile, &
    t2bulkdensity,t2moisture,t2ss,tcrowndiameter,tunique_species
  use io_variables, only : verbose,workdir,filesep
  use constant_variables, only : tolerance
  implicit none
  
  ! Local Variables
  integer :: i,j,it,tindex,itree,ierror
  real :: nsub,rnumx,rnumy,newx,newy
  real :: treeid,xtest,ytest
  real :: x_loc_min,x_loc_max,y_loc_min,y_loc_max,x_per,y_per
  real :: dataleft,dataright,databottom,datatop
  character(len=50) :: treelistformat
  integer,allocatable :: numarray(:)
  integer,allocatable :: rounddown(:),ntreesold(:)
  real,allocatable :: temp_array(:)
  
  ! Variables for finding height to max crown radius
  real:: beta_a, beta_b, beta_c, beta_norm, z, rad, norm_height
  integer:: k
  integer,dimension(2):: max_ind
  real,parameter,dimension(11) :: z_range = &
    (/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
  real, dimension(2,11):: z_rad
  
  ! Executable Code
  ! Find and trim locations and dimensions from FF treelist
  if(itrees.eq.7)then ! Find data and trim data locations
    open(48,file=TRIM(TRIM(workdir)//filesep)//treefile)
    !fast fuels list x/y locations as coordinates, find min/max of x/y_loc_max
    !using x/y extremes, shift xpos and ypos
    x_loc_max = 0
    y_loc_max = 0
    x_loc_min = 0
    y_loc_min = 0
    
    allocate(temp_array(19))
    read(48,*) !read 1st line and throw away, has column headers
    do
      read(48,*,iostat=ierror) temp_array(:)
      if(ierror<0)then
        close(48)
        exit
      endif
      if (i.eq.1) then  
        x_loc_max = temp_array(18) !initilize min/max
        x_loc_min = temp_array(18) !initilize min/max
        y_loc_max = temp_array(19) !initilize min/max
        y_loc_min = temp_array(19) !initilize min/max
      endif
      x_loc_max = max(temp_array(18), x_loc_max) !max x for transformation
      x_loc_min = min(temp_array(18), x_loc_min) !min x for transformation
      y_loc_max = max(temp_array(19), y_loc_max) !max y for transformation
      y_loc_min = min(temp_array(19), y_loc_min) !min y for transformation
    enddo
    print*, 'x min: ',x_loc_min,' x max: ',x_loc_max
    print*, 'y min: ',y_loc_min,' y max: ',y_loc_max
    
    x_per = nint(sqrt( (x_loc_min - x_loc_max)**2 ))
    y_per = nint(sqrt( (y_loc_min - y_loc_max)**2 ))
    
    !cut dataset to negate edge effects, use ndatax, ndatay
    !if larger dataset specified that what exists, revert to 
    !the largest part of data
    ndatax = x_per
    print*,'New data x length: ', ndatax
    ndatay = y_per
    print*,'New data y length: ', ndatay
    datalocx = 0
    datalocy = 0
  endif

  !---Determine how many dataset subdomains fit within your main domain
  nsub = (nx*dx*ny*dy)/(ndatax*ndatay)
  print*,'Number of subdomains = ',nsub
  allocate(ntreesold(ntspecies))
  ntreesold = ntrees
  ntrees=nint(ntrees*nsub)
  print*,'Additional trees needed to fill domain = ',ntrees-ntreesold
 
  !---Allocate and fill arrays 
  allocate(tlocation(ntspecies,maxval(ntrees),2))     ! Tree cartesian coordinates [m,m]
  allocate(theight(maxval(ntrees),ntspecies))         ! Tree heights [m]
  allocate(tcrownbotheight(maxval(ntrees),ntspecies)) ! Height to live crown [m]
  allocate(tcrowndiameter(maxval(ntrees),ntspecies))  ! Crown diameter [m]
  allocate(tcrownmaxheight(maxval(ntrees),ntspecies)) ! Height to max crown diameter [m]
  
  allocate(t2bulkdensity(maxval(ntrees),tfuelbins,ntspecies)) ! Crown fuel bulk density [kg/m3]
  allocate(t2moisture(maxval(ntrees),tfuelbins,ntspecies))    ! Crown fuel moisture content [fraction]
  allocate(t2ss(maxval(ntrees),tfuelbins,ntspecies))          ! Crown fuel sizescale [m]
  allocate(numarray(ntspecies)); numarray(:)=0
  
  !---Fill arrays wit tree list data
  open (48,file=TRIM(TRIM(workdir)//filesep)//treefile, &
    form='formatted',status='old')
  if(itrees.eq.7)then ! Fast fuels tree list
    read(48,*)
    do
      read(48,*,iostat=ierror) temp_array(:)
      if(ierror<0)then
        close(48)
        exit
      endif
      tindex=findloc(tunique_species,temp_array(5),dim=1)
      numarray(tindex) = numarray(tindex)+1
      ! Tree location
      tlocation(tindex,numarray(tindex),1)=abs(temp_array(18)-x_loc_min) !x loc
      tlocation(tindex,numarray(tindex),2)=abs(temp_array(19)-y_loc_min) !y loc
      ! Height, HTLC, crown diameter
      theight(numarray(tindex),tindex)         = temp_array(3)!/3.281 !convert ft to meters
      tcrownbotheight(numarray(tindex),tindex) = temp_array(10)!/3.281 !convert ft to meters
      tcrowndiameter(numarray(tindex),tindex)  = temp_array(17)*2 !crown rad * 2
      ! Find crown max height
      beta_a = temp_array(11)
      beta_b = temp_array(12)
      beta_c = temp_array(13)
      beta_norm = temp_array(14)
      do k=1,11
        z = z_range(k)
        rad = (beta_c*(z**(beta_a-1)*(1-z)**(beta_b-1)))/beta_norm
        z_rad(1,k) = z
        z_rad(2,k) = rad
      enddo
      max_ind = maxloc(z_rad, dim=2)
      norm_height = z_rad(1,max_ind(2))
      tcrownmaxheight(numarray(tindex),tindex) = &
        ((norm_height*temp_array(9))+temp_array(10))!/3.281 !convert ft to meters
      ! Bulk density, moisture, size scale
      do j=1,tfuelbins
        t2bulkdensity(numarray(tindex),j,tindex) = &
          temp_array(15)/temp_array(16) ! weight (kg) / volume (m^3)
        t2moisture(numarray(tindex),j,tindex) = 1.0 ! setting canopy moisture to 100% ; not avalible in FF
        t2ss(numarray(tindex),j,tindex) = min(0.002,2./temp_array(8)) ! suppose infinite cylinder, (2*pi*r*L)/(pi*r^2*L) = sav, solve for r
      enddo
    enddo
  else ! Forest Service tree list
    allocate(temp_array(7+3*tfuelbins))
    do
      read(48,*,iostat=ierror) temp_array(:)
      if(ierror<0)then
        close(48)
        exit
      endif
      tindex=findloc(tunique_species,temp_array(1),dim=1)
      numarray(tindex) = numarray(tindex)+1
      if (itrees.eq.2) then
        tlocation(tindex,numarray(tindex),1) = &
          temp_array(2)+datalocx
        tlocation(tindex,numarray(tindex),2) = &
          temp_array(3)+datalocy
      else
        call random_number(xtest)
        xtest = ndatax*xtest
        call random_number(ytest)
        ytest = ndatay*ytest
    
        tlocation(tindex,numarray(tindex),1) = xtest+datalocx
        tlocation(tindex,numarray(tindex),2) = ytest+datalocy
      endif
      theight(numarray(tindex),tindex) = temp_array(4)
      tcrownbotheight(numarray(tindex),tindex) = temp_array(5)
      tcrowndiameter(numarray(tindex),tindex) = temp_array(6)
      tcrownmaxheight(numarray(tindex),tindex) = temp_array(7)
      do j=1,tfuelbins
        t2bulkdensity(numarray(tindex),j,tindex) = &
          temp_array(8+3*(j-1))
        t2moisture(numarray(tindex),j,tindex) = &
          temp_array(9+3*(j-1))
        t2ss(numarray(tindex),j,tindex) = temp_array(10+3*(j-1))
      enddo
    enddo
  endif
  deallocate(numarray,temp_array)
  
  !---Finish filling arrays with replicas of treelist data
  ! For each of the trees in the dataset which are placed above, we copy the
  ! information for each of them into the arrays and choose a new location;
  ! if the initial random location is within the chosen dataset area, we
  ! choose another one
  if(nsub.gt.1) then
    ! First we need to find the area for where the dataset will live
    dataleft = datalocx
    dataright = (datalocx + ndatax)
    databottom = datalocy
    datatop = (datalocy + ndatay)
  
    print*,'dataset lives at these coordinates: ', &
      dataleft,dataright,databottom,datatop

    do i=1,ntspecies
      do j=1,ntrees(i)-ntreesold(i) 
        ! Sample from tree list for new tree
        call random_number(treeid)
        tindex=floor(treeid*(ntrees(i)+1))

        ! Find location for new tree
        it=j+ntreesold(i)
        newx = tlocation(i,tindex,1)
        newy = tlocation(i,tindex,2)
        do while (newx.ge.dataleft.and.newx.le.dataright.and. &
          newy.ge.databottom.and.newy.le.datatop)
          call random_number(rnumx)
          xtest = rnumx*nx*dx
          call random_number(rnumy)
          ytest = rnumy*ny*dy
          if(any(sqrt((xtest-tlocation(:,1:it-1,1))**2+ &
            (ytest-tlocation(:,1:it-1,2))**2).lt.0.1))then
            cycle
          else
            newx=xtest
            newy=ytest
          endif
        enddo
        tlocation(i,it,1) = newx
        tlocation(i,it,2) = newy
        theight(it,i) = theight(tindex,i)
        tcrownbotheight(it,i) = tcrownbotheight(tindex,i)
        tcrowndiameter(it,i) = tcrowndiameter(tindex,i)
        tcrownmaxheight(it,i) = tcrownmaxheight(tindex,i)
        t2bulkdensity(it,:,i) = t2bulkdensity(tindex,:,i)
        t2moisture(it,:,i) = t2moisture(tindex,:,i)
        t2ss(it,:,i) = t2ss(tindex,:,i)
      enddo
    enddo
  endif
  
  if (verbose.eq.1) then 
    call define_newtreefile
    open (222,file=TRIM(newtreefile)//'_treelist.txt', &
      form='formatted',status='unknown')
    treelistformat = '(I2.1,6F12.5,F10.6,F10.4,F10.5)'
    do j=1,ntspecies
      do i=1,ntrees(j)
        !species, x, y, ht, htlc, cd , htmcd, cbd, fmc, ss
        write(222,treelistformat) j,&
          tlocation(j,i,1),&
          tlocation(j,i,2),&
          theight(i,j),&
          tcrownbotheight(i,j),&
          tcrowndiameter(i,j),&
          tcrownmaxheight(i,j),&
          t2bulkdensity(i,1,j),&
          t2moisture(i,1,j),&
          t2ss(i,1,j)
      enddo
    enddo
    close(222)
  endif
end subroutine treelist_readin

!-----------------------------------------------------------------------
! This subroutine will find the number of species in any treelist
!-----------------------------------------------------------------------
subroutine find_numspecies
  use fuel_variables, only : treefile,itrees,iFIAspecies,ntspecies, &
    ntrees,tunique_species
  use io_variables, only : workdir,filesep
  implicit none
  
  !Local Variables
  integer i,numtrees,min_val_sp,max_val_sp,ierror
  integer,allocatable :: unique_species(:)
  integer,allocatable :: temp_array(:)
  real,dimension(19) :: read_array
  
  ! Executable Code
  ! Count number of trees in file
  numtrees = 0
  open (2,file=TRIM(TRIM(workdir)//filesep)//treefile,status= 'old')
  do
    read (2,*,iostat=ierror)
    if(ierror<0)then
      rewind(2)
      exit
    else
      numtrees=numtrees+1
    endif
  enddo
  
  ! Find the number of tree species in file
  if(itrees.eq.7)then ! Read FastFuels Files
    numtrees=numtrees-1
    allocate(temp_array(numtrees))
    read(2,*) !read 1st line and throw away, has column headers
    do i=1,numtrees
      read(2,*) read_array(:)
      if (iFIAspecies.eq.1) then
        temp_array(i)=int(read_array(1))
      else
        temp_array(i)=int(read_array(5)) !take from sp_grp, 5th pos
      endif
    enddo
  else
    allocate(temp_array(numtrees))
    do i=1,numtrees
      read(2,*) temp_array(i)
    enddo
  endif
  close(2)
  
  ! Organize species
  allocate(unique_species(maxval(temp_array)))
  min_val_sp = minval(temp_array)-1
  max_val_sp = maxval(temp_array)
  
  i=0
  do while (min_val_sp<max_val_sp)
     i=i+1
     min_val_sp=minval(temp_array,mask=temp_array>min_val_sp)
     unique_species(i)=min_val_sp
  enddo
  ntspecies=i
  allocate(tunique_species(i),source=unique_species(1:i))
  
  ! Count number of trees in each species
  allocate(ntrees(ntspecies))
  ntrees=0
  do i=1,numtrees
    ntrees(findloc(tunique_species,temp_array(i)))= &
      ntrees(findloc(tunique_species,temp_array(i)))+1
  enddo
  deallocate(temp_array,unique_species)
  
  print*,'Number of Set Species = ',ntspecies
  print*,'Groups = ',tunique_species
  print*,'Number of Trees in each Species = ',ntrees
  
end subroutine find_numspecies
