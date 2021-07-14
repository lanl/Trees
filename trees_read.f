
      subroutine treesGeneral_readin
      !-----------------------------------------------------------------
      ! treesGeneral_readin is a function which reads in a general trees
      ! file for use in FIRETEC or QUIC-Fire
      !-----------------------------------------------------------------
      use grid_variables
      use baseline_variables
      
      implicit none
      integer i
        
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
      
      end subroutine treelist_readin
      
      subroutine json_readin
      !-----------------------------------------------------------------
      ! json_readin is a function which reads in .json files for use in 
      ! FIRETEC or QUIC-Fire
      !-----------------------------------------------------------------
      use grid_variables
      use baseline_variables
      use json_module
      implicit none
      type(json_file) :: json
      type(json_core) :: jCore
      character(LEN=20):: istr1,istr2
      logical :: isFound
      integer :: i,j
      integer :: nplots,ttemp,temp
      real    :: rtemp,xtemp,ytemp

      ! Initialize the json_file object
      call json%initialize()
      
      ! Load the file
      call json%load_file('treelist.json'); if(json%failed()) stop

      ! Determine number of trees
      ntspecies=1
      allocate(ntrees(ntspecies))
      call json%get_core(jCore)
      call json%info('plots', n_children=nplots)
      ntrees(1)=0
      do i=1,nplots
        write(istr1,'(I5)') i
        call json%info('plots('//istr1//').trees', n_children=temp)
        ntrees(1)=ntrees(1)+temp
      enddo
      print*,'JSON file contains information for',ntrees(1),'trees'

      ! Allocate needed arrays
      allocate(tlocation(ntspecies,maxval(ntrees),2)) ! Tree cartesian coordinates [m,m]
      allocate(theight(maxval(ntrees),ntspecies)) ! Tree heights [m]
      allocate(tcrownbotheight(maxval(ntrees),ntspecies)) ! Height to live crown [m]
      allocate(tcrowndiameter(maxval(ntrees),ntspecies)) ! Crown diameter [m]
      allocate(tcrownmaxheight(maxval(ntrees),ntspecies)) ! Height to max crown diameter [m]
      allocate(t2bulkdensity(maxval(ntrees),tfuelbins,ntspecies)) ! Crown fuel bulk density [kg/m3]
      allocate(t2moisture(maxval(ntrees),tfuelbins,ntspecies)) ! Crown fuel moisture content [fraction]
      allocate(t2ss(maxval(ntrees),tfuelbins,ntspecies)) ! Crown fuel size scale [m]

      ! Read in data
      temp=0
      do i=1,nplots
        write(istr1,'(I5)') i
        call json%get('plots('//istr1//').x',xtemp,isFound)
        call json%get('plots('//istr1//').y',ytemp,isFound)
        call json%info('plots('//istr1//').trees', n_children=ttemp)
        do j=1,ttemp
          temp=temp+1
          write(istr1,'(I5)') i
          write(istr2,'(I5)') j
          call json%get('plots('//istr1//').trees('//istr2//').x',rtemp,isFound)
          tlocation(1,temp,1)=xtemp+rtemp
          call json%get('plots('//istr1//').trees('//istr2//').y',rtemp,isFound)
          tlocation(1,temp,2)=ytemp+rtemp
          call json%get('plots('//istr1//').trees('//istr2//').height',theight(temp,1),isFound)
          call json%get('plots('//istr1//').trees('//istr2//').crownBaseHeight',tcrownbotheight(temp,1),isFound)
          call json%get('plots('//istr1//').trees('//istr2//').crownRadius',rtemp,isFound)
          tcrowndiameter(temp,1)=2.*rtemp
          call json%get('plots('//istr1//').trees('//istr2//').crownHeight',rtemp,isFound)
          tcrownmaxheight(temp,1)=tcrownbotheight(temp,1)+rtemp/3.
          call json%get('plots('//istr1//').trees('//istr2//').crownBulkDensity',t2bulkdensity(temp,1,1),isFound)
          t2moisture(temp,1,1)=t2bulkdensity(temp,1,1)
          call json%get('plots('//istr1//').trees('//istr2//').inverseSav',rtemp,isFound)
          t2ss(temp,1,1)=2./rtemp
        enddo
      enddo
      xtemp=minval(tlocation(1,:,1))
      ytemp=minval(tlocation(1,:,2))
      do i=1,ntrees(1)
        tlocation(1,i,1)=tlocation(1,i,1)-xtemp+5
        tlocation(1,i,2)=tlocation(1,i,2)-ytemp+5
      enddo
      print*,"Size of fuel domain:"
      print*,"x(min)",minval(tlocation(1,:,1)),"x(max)",maxval(tlocation(1,:,1))
      print*,"y(min)",minval(tlocation(1,:,2)),"y(max)",maxval(tlocation(1,:,2))
      
      end subroutine json_readin
