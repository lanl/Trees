!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! treatments contains functions which prescribe different treatments
! performed on the baseline forest
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine treatment
      !-----------------------------------------------------------------
      ! treatment is a function which calls the appropriate treatment
      ! functions to fill the rhof, sizescale, moisture, and
      ! fueldepth arrays.
      !-----------------------------------------------------------------
      use grid_variables
      use treatment_variables

      implicit none
      
      if (itreatment.eq.1) call slash_layer
      if (itreatment.eq.2) call slash_piles
      if (itreatment.eq.3) call clearing

      end subroutine treatment

      subroutine slash_layer
!----------------------------------------------------------------------
! slash_layer is a function takes all of the canopy fine fuels and lays
! them on the ground in a layer with thickness defined in the variables
! file.
!----------------------------------------------------------------------
      use grid_variables
      use treatment_variables
      use baseline_variables

      implicit none
      
      integer i,j,k,ift
      real Cremove
      real,allocatable:: srho(:)
      real,allocatable:: total_mass(:),total_ss(:)
      real actual_mass,target_mass

      allocate(total_mass(nfuel))
      allocate(total_ss(nfuel))
      allocate(srho(nfuel))
      
      print*,'Creating Slash Layer'
      ! Sum together all the mass throughout the entire domain
      target_mass = 0
      do ift=1,nfuel
        do i=1,nx
          do j=1,ny
            do k=1,zmax
              target_mass = target_mass+rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
            enddo
          enddo
        enddo
      enddo
      total_mass(:) = 0.
      total_ss(:)   = 0.
      do ift=1,nfuel
        do i=int(sdnx(1)/dx+1),int(sdnx(2)/dx+1)
          do j=int(sdny(1)/dy+1),int(sdny(2)/dy+1)
            do k=1,zmax
              Cremove = 0.
              if (istem.eq.1)then
                if (MOD(ift-ngrass,ntreefueltypes).eq.tfuelbins+1) Cremove = 0.8
                if (MOD(ift-ngrass,ntreefueltypes).eq.tfuelbins+2) Cremove = 0.4
              endif
              total_mass(ift) = total_mass(ift)+(1-Cremove)*rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
              total_ss(ift)  = total_ss(ift)+sizescale(ift,i,j,k)*(1-Cremove)*rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
              target_mass = target_mass-Cremove*rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
            enddo
          enddo
        enddo

        ! Compute mass density of slash layer give
        srho(ift) = total_mass(ift)/((sdnx(2)-sdnx(1))*(sdny(2)-sdny(1))*sdepth)
      enddo
      print*,'Mass density of slash layer',sum(srho)

      ! Distribute that mass into the array
      do ift=1,nfuel
        do i=int(sdnx(1)/dx+1),int(sdnx(2)/dx+1)
          do j=int(sdny(1)/dy+1),int(sdny(2)/dy+1)
            fueldepth(ift,i,j,1)= sdepth
            do k=1,zmax
              sizescale(ift,i,j,k) = total_ss(ift)/total_mass(ift)
              if(zheight(i,j,k+1).lt.sdepth) then
                rhof(ift,i,j,k) = srho(ift)
                moist(ift,i,j,k)= smoist
              else if(zheight(i,j,k).lt.sdepth) then
                rhof(ift,i,j,k) = srho(ift)*(sdepth-zheight(i,j,k))/(zheight(i,j,k+1)-zheight(i,j,k))
                moist(ift,i,j,k)= smoist
              else
                rhof(ift,i,j,k) = 0.
                moist(ift,i,j,k)= 0.
              endif
            enddo
          enddo
        enddo
      enddo
      
      ! Print out the target and actual fuel masses for comparisons sake
      print*,'Slash layer target fuel mass:',target_mass
      actual_mass = 0
      do ift=1,nfuel
        do i=1,nx
          do j=1,ny
            do k=1,zmax
              actual_mass = actual_mass+rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
            enddo
          enddo
        enddo
      enddo
      print*,'Slash layer actual fuel mass:',actual_mass
      print*,'Slash layer error:',actual_mass/target_mass*100,'%'

      end subroutine slash_layer

      subroutine slash_piles
!----------------------------------------------------------------------
! slash_pile is a function takes all of the canopy fine fuels and
! stacks them in piles with dimensions defined in the variables
! file.
!----------------------------------------------------------------------
      use constant_variables
      use grid_variables
      use treatment_variables
      use baseline_variables

      implicit none
      
      integer ift,i,ii,iii,j,jj,jjj,k,kk,kkk
      integer snx,sny,sxnum,synum
      real Cremove
      real snumber
      real xtest,ytest
      integer xtestbot,xtesttop,ytestbot,ytesttop,ztesttop
      real aellip,rhoftemp
      integer cellcount
      real spmass,spheight
      real target_mass,actual_mass
      real,allocatable:: total_mass(:),total_ss(:)
      real,external:: paraboloid
      
      allocate(total_mass(nfuel))
      allocate(total_ss(nfuel))

      print*,'Creating Slashpiles'
      ! Sum together all the mass throughout the entire domain minus the ground layer
      target_mass = 0
      do ift = 1,nfuel
        if (ift.gt.ngrass) then
          do i=1,nx
            do j=1,ny
              do k=1,zmax
                target_mass = target_mass+rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
              enddo
            enddo
          enddo
          total_ss(ift)   = 0.
          total_mass(ift) = 0.
          do i=int(sdnx(1)/dx+1),int(sdnx(2)/dx+1)
            do j=int(sdny(1)/dy+1),int(sdny(2)/dy+1)
              do k=1,zmax
                Cremove = 0.
                if (istem.eq.1)then
                  if (MOD(ift-ngrass,ntreefueltypes).eq.tfuelbins+1) Cremove = 0.8
                  if (MOD(ift-ngrass,ntreefueltypes).eq.tfuelbins+2) Cremove = 0.4
                endif
                total_mass(ift)= total_mass(ift)+(1-Cremove)*rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
                total_ss(ift)  = total_ss(ift)+sizescale(ift,i,j,k)*(1-Cremove)*rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
                target_mass = target_mass-Cremove*rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
                rhof(ift,i,j,k) = 0.
                moist(ift,i,j,k)= 0.
                sizescale(ift,i,j,k) = 0.
              enddo
            enddo
          enddo
        endif
      enddo

      ! Compute the number and location of slash piles within the clearing
      snx = sdnx(2)-sdnx(1)
      sny = sdny(2)-sdny(1)
      snumber = sum(total_mass)/(PI/4.*sdiameter**2.*sheight*sprho)
      print*,'Total number of slash piles:',ceiling(snumber)
      sxnum   = floor(snx/(snx*sny/ceiling(snumber))**0.5)
      synum   = ceiling(snumber/sxnum)
      ! Resolve the full slash piles
      do i=1,floor(snumber)
        xtest = sdnx(1)+(1+mod(i-1,sxnum))*snx/(sxnum+1)
        xtestbot = floor((xtest-sdiameter/2.)/dx)
        xtesttop = floor((xtest+sdiameter/2.)/dx)
        ytest = sdny(1)+(1+((i-1)/sxnum))*sny/(synum+1)
        ytestbot = floor((ytest-sdiameter/2.)/dx)
        ytesttop = floor((ytest+sdiameter/2.)/dx)
        do k=1,zmax
          if(zheight(nint(xtest/dx),nint(ytest/dy),k).gt.sheight) then
            ztesttop = k
            exit
          endif
        enddo
        
        aellip = -sdiameter**2./(2.*sheight)
        do ii=xtestbot,xtesttop
          do jj=ytestbot,ytesttop
            do kk=1,ztesttop
              ! Determine how many of subcells of a cell are within the paraboloid, the fraction of the subcells is equal to the fraction of the cell within the paraboloid
              cellcount = 0
              do iii=1,10
                do jjj=1,10
                  do kkk=1,10
                    if(paraboloid(aellip,(ii+(2.*iii-1.)/20.)*dx,xtest,(jj+(2.*jjj-1)/20.)*dy,ytest,sheight).gt.
     +                zheight(ii,jj,kk)+(2.*kkk-1.)/20.*(zheight(ii,jj,kk+1)-zheight(ii,jj,kk)))  cellcount = cellcount+1
                  enddo
                enddo
              enddo
              ! Fill in the 3D arrays for full slashpiles
              do ift=1,nfuel
                rhoftemp = sprho*cellcount/1000.*total_mass(ift)/sum(total_mass)
                if(rhoftemp.gt.0) then
                  moist(ift,ii,jj,kk) = (rhoftemp*smoist+rhof(ift,ii,jj,kk)*moist(ift,ii,jj,kk))/(rhoftemp+rhof(ift,ii,jj,kk))
                  rhof(ift,ii,jj,kk)  = rhof(ift,ii,jj,kk) + rhoftemp
                  sizescale(ift,ii,jj,kk) = total_ss(ift)/total_mass(ift)
                endif
              enddo
            enddo
          enddo
        enddo
      enddo

      ! Resolve the last partial pile (same diameter different height)
      if (ceiling(snumber).eq.1) i=0
      ! Hack in the case that there's only 1 pile
      spmass  = sum(total_mass)-floor(snumber)*(PI/4.*sdiameter**2.*sheight*sprho)
      spheight= 3.*spmass/(sprho*PI*sdiameter**2.)
      xtest = sdnx(1)+(1+mod(i-1,sxnum))*snx/(sxnum+1)
      xtestbot = floor((xtest-sdiameter/2.)/dx)
      xtesttop = floor((xtest+sdiameter/2.)/dx)
      ytest = sdny(1)+(1+(i-1)/sxnum)*sny/(synum+1)
      ytestbot = floor((ytest-sdiameter/2.)/dx)
      ytesttop = floor((ytest+sdiameter/2.)/dx)
      do k=1,zmax
        if(zheight(nint(xtest/dx),nint(ytest/dy),k+1).gt.sheight) then
          ztesttop = k
          exit
        endif
      enddo
     
      aellip = -sdiameter**2./(2.*spheight)
      do ii=xtestbot,xtesttop
        do jj=ytestbot,ytesttop
          do kk=1,ztesttop
            ! Determine how many of subcells of a cell are within the paraboloid, the fraction of the subcells is equal to the fraction of the cell within the paraboloid
            cellcount = 0
            do iii=1,10
              do jjj=1,10
                do kkk=1,10
                  if(paraboloid(aellip,(ii+(2.*iii-1.)/20.)*dx,xtest,(jj+(2.*jjj-1)/20.)*dy,ytest,spheight).gt.
     +              zheight(ii,jj,kk)+(2.*kkk-1.)/20.*(zheight(ii,jj,kk+1)-zheight(ii,jj,kk)))  cellcount = cellcount+1
                enddo
              enddo
            enddo

            ! Fill in the 3D arrays for partial slashpile
            do ift=1,nfuel
              rhoftemp = sprho*cellcount/1000.*total_mass(ift)/sum(total_mass)
              if(rhoftemp.gt.0) then
                moist(ift,ii,jj,kk) = (rhoftemp*smoist+rhof(ift,ii,jj,kk)*moist(ift,ii,jj,kk))/(rhoftemp+rhof(ift,ii,jj,kk))
                rhof(ift,ii,jj,kk)  = rhof(ift,ii,jj,kk) + rhoftemp
                sizescale(ift,ii,jj,kk) = total_ss(ift)/total_mass(ift)
              endif
            enddo
          enddo
        enddo
      enddo

      ! Print out the target and actual fuel masses for comparisons sake
      print*,'Slashpile target fuel mass:',target_mass
      actual_mass = 0
      do ift=1,nfuel
        do i=1,nx
          do j=1,ny
            do k=1,zmax
              actual_mass = actual_mass+rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
            enddo
          enddo
        enddo
      enddo
      print*,'Slashpile actual fuel mass:',actual_mass
      print*,'Slashpile error:',actual_mass/target_mass*100,'%'
    
      end subroutine slash_piles
      
      subroutine clearing
!----------------------------------------------------------------------
! clearing is a function removes all fuels in the treatment area.
!----------------------------------------------------------------------
      use grid_variables
      use treatment_variables

      implicit none
      
      integer i,j,k,ift
      
      ! Clear fuels out of the treatment area
      do ift=1,nfuel
        do i=int(sdnx(1)/dx+1),int(sdnx(2)/dx+1)
          do j=int(sdny(1)/dy+1),int(sdny(2)/dy+1)
            fueldepth(ift,i,j,1) = zheight(i,j,k+1)
            do k=1,zmax
              sizescale(ift,i,j,k) = 0
              rhof(ift,i,j,k) = 0
              moist(ift,i,j,k)= 0
            enddo
          enddo
        enddo
      enddo

      end subroutine clearing
