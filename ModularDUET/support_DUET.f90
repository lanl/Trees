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
!-----------------------------------------------------------------

! This is a list of subroutines for the DUET program.
! Author:  Jenna Sjunneson McDanold 7/2024


module support
  implicit none
  contains
    subroutine decay
        use DUETio

        integer :: s,i,j,yt,yr,spy

        print*,'Decay and compaction occurring...'
        print*,'max and min of lrho: ',maxval(litter%lrho),minval(litter%lrho)
        ! Decay and compact the litter
        do s=1,domain%ns
        
          do i=1,domain%nx
            do j=1,domain%ny
              yt=1
              do yr=1,domain%ysb
                do spy=1,domain%spy
                  litter%lrho(s,i,j,yt) = litter%lrho(s,i,j,yt)*exp(-species%decy(s)*(yt))
                yt=yt+1
                enddo
              enddo
            enddo
          enddo
        enddo
        print*,'max and min of lrho after decay: ',maxval(litter%lrho),minval(litter%lrho)
                
    end subroutine decay

    !---------------------------------------------------------------------------------------------------
    
    subroutine diffuse

      use DUETio

      integer :: i,j,s,yt,ct
      integer :: iplus,iminus,jplus,jminus
      real :: a,deltat,mval !,densitythresh
      real,allocatable :: tmp(:,:,:,:), TEMP(:,:,:,:)
      
      print*,'Sum of litter before diffusion: ',sum(litter%lrho)
      print*,'Minimum possible thresh value:',sum(litter%lrho)/(domain%nx*domain%ny)
      print*,'Density threshold is:',duetvars%densitythresh
      
      if((sum(litter%lrho)/(domain%nx*domain%ny)).ge.duetvars%densitythresh) then
          print*,'DANGER WILL ROBINSON!'
          print*,'DENSITYTHRESH IS TOO LOW!!!  MUST BE GREATER THAN ',sum(outarray%lrho)/(domain%nx*domain%ny)
          print*,'ADJUSTING...'
          duetvars%densitythresh = (1.1)*sum(outarray%lrho)/(domain%nx*domain%ny)
          print*,'Densitythresh reset to ',duetvars%densitythresh
      endif
              
      a = (domain%dx**2)/4
      deltat = (domain%dy**2)/4
      
      allocate(tmp(domain%ns,domain%nx,domain%ny,domain%nt))
      allocate(TEMP(domain%ns,domain%nx,domain%ny,domain%nt))
      tmp = 0.0
      TEMP = 0.0
      
      ! Diffusion from each cell to the surrounding cells provided the density is less
      ! than that within the center cell
      
      tmp = litter%lrho
      
      ct = 0
      if(maxval(sum(tmp,DIM=1)).gt.duetvars%densitythresh) then
        do while (maxval(sum(tmp,DIM=1)).gt.duetvars%densitythresh)
          mval = maxval(sum(tmp,DIM=1))
          do j=1,domain%ny
            do i=1,domain%nx
              do s=1,domain%ns
                do yt=1,domain%nt
                  iplus = i+1
                  jplus = j+1
                  iminus = i-1
                  jminus = j-1
      
                  if (i.eq.1) iminus = domain%nx
                  if (i.eq.domain%nx) iplus = 1
                  if (j.eq.1) jminus = domain%ny
                  if (j.eq.domain%ny) jplus = 1
      
                      ! 2D diffusion equation
                    TEMP(s,i,j,yt) = a*deltat*((tmp(s,iplus,j,yt)-2*tmp(s,i,j,yt)+tmp(s,iminus,j,yt))/domain%dx**2 &
                      + (tmp(s,i,jplus,yt)-2*tmp(s,i,j,yt)+tmp(s,i,jminus,yt))/domain%dy**2) &
                      + tmp(s,i,j,yt)
                enddo
              enddo
            enddo
          enddo
          
          tmp = TEMP
          ct = ct+1
          if(mval-maxval(sum(TEMP,DIM=1)).lt.0.0001) then
            print*,'Reached steady state with litter diffusion.'
            print*,'Max of litter =',mval
            exit
          endif
        enddo
      else
        TEMP = tmp
      endif
      print*,'Diffusion complete'
      print*,'Diffusion ran for',ct,'rounds'
      print*,'Sum after diffusion = ',sum(TEMP)

      litter%lrho = TEMP        

    end subroutine diffuse

    !---------------------------------------------------------------------------------------------------
    
    subroutine depthMoistSs

      use DUETio

      integer :: s,i,j,yt,yr,spy

      do s=1,domain%ns
        do i=1,domain%nx
          do j=1,domain%ny
            yt=1
            do yr=1,domain%ysb
              do spy=1,domain%spy
                if (litter%lrho(s,i,j,yt).gt.0) then 
                  litter%lafd(s,i,j,yt) = min((litter%lrho(s,i,j,yt)/(domain%dx*domain%dy) &
                    *exp(species%comp(s)*(yt))),inarray%zheight(i,j,2)-inarray%zheight(i,j,1))
                  litter%lh20(s,i,j,yt) = species%fh20(s)/100*exp(-2.0*yt)
                  if (litter%lh20(s,i,j,yt).lt.duetvars%relhum) litter%lh20(s,i,j,yt) = duetvars%relhum
                  if (litter%lrho(s,i,j,yt).gt.0) litter%lsss(s,i,j,yt) = species%ssss(s)
                endif 
                yt=yt+1
              enddo
            enddo
            outarray%lrho(s,i,j) = sum(litter%lrho(s,i,j,:))
            outarray%lafd(s,i,j) = maxval(litter%lafd(s,i,j,:)/domain%ns)     
            outarray%lh20(s,i,j) = maxval(litter%lh20(s,i,j,:)/domain%ns)
            outarray%lsss(s,i,j) = maxval(litter%lsss(s,i,j,:))     
          enddo
        enddo
      enddo

    end subroutine depthMoistSs

    !---------------------------------------------------------------------------------------------------
    
    subroutine growGrass

        use DUETio

        integer :: i,j,k,s,zmax,ct,y,yt
        real :: g 
        real :: litterfactor,shadefactor,rhocolumn
        real,allocatable :: grhof(:,:,:,:)

        allocate(grhof(domain%ng,domain%nx,domain%ny,domain%nt))

        zmax = domain%nz-1

        print*,'Grass growing...'

        ! Add grass
        do s=1,domain%ng
          do i=1,domain%nx
            do j=1,domain%ny
              rhocolumn = 0
              do k=1,zmax
                rhocolumn = rhocolumn+sum(inarray%trhof(:,i,j,k))* &
                  (inarray%zheight(i,j,k+1)-inarray%zheight(i,j,k))/inarray%zheight(i,j,zmax+1)
              enddo
              shadeFactor = exp(-duetvars%grassconstant*rhocolumn/0.6)
              litterFactor=exp(-duetvars%litterconstant*sum(outarray%lrho(:,i,j))/0.6)
              ct=1
              do y=1,domain%ysb
                  do yt=1,domain%spy
                      ! For every cell that includes grass we will determine the height of the
                      ! grass, the amount of shade from the canopy, the the amount of cover from
                      ! the litter, and then build the grass density   

                      call random_number(g)
                      if(g.lt.0.75) g = g+0.75
                      if(g.gt.1) g = 1
                      grhof(s,i,j,ct)=g*grasses%grho(s)*shadeFactor*litterFactor &
                        *exp(-grasses%decy(s)*(ct-1))
                      ct=ct+1
                  enddo
              enddo
            enddo
          enddo
        enddo      
        print*,'Max, min of grass = ',maxval(grhof(1,:,:,:)),minval(grhof(1,:,:,:))
        print*,'Total grass = ',sum(grhof(1,:,:,:))

        do s=1,domain%ng
          do i=1,domain%nx
            do j=1,domain%ny
              outarray%srho(s,i,j) = sum(grhof(s,i,j,:))
              outarray%frho(s,i,j,1) = sum(grhof(s,i,j,:))
              if(outarray%srho(s,i,j).gt.0) then
                outarray%safd(s,i,j) = grasses%dept(s)*outarray%srho(s,i,j)
                outarray%ssss(s,i,j) = grasses%ssss(s)
                outarray%sh20(s,i,j) = grasses%fh20(s)*outarray%srho(s,i,j)
                outarray%fafd(s,i,j) = grasses%dept(s)*outarray%srho(s,i,j)
                outarray%fsss(s,i,j,1) = grasses%ssss(s)
                outarray%fh20(s,i,j,1) = grasses%fh20(s)*outarray%srho(s,i,j)
              endif
            enddo
          enddo
        enddo
        print*,'Max, min of grass moisture = ',maxval(outarray%sh20),minval(outarray%sh20)
        deallocate(grhof)

      
    end subroutine growGrass


    !subroutine grassmoistureoverride
    !  use duetio
!
    !  if (duetvars%gmoistoverride.ne.0) then
    !    print*,'gmoistoverride = ',duetvars%gmoistoverride
    !    print*,'max value of moisture = ',maxval(outarray%lh20(:,:,:)+moist(ngrass+fueltotal:,:,:,1))
    !    outarray%lh20(:,:,:) = duetvars%gmoistoverride/maxval(outarray%lh20(:,:,:)+moist(ngrass+fueltotal:,:,:)) * outarray%lh20(:,:,:)
    !    moist(ngrass+fueltotal:,:,:,1) = duetvars%gmoistoverride/maxval(outarray%lh20(:,:,:)+ &
    !      moist(ngrass+fueltotal:,:,:,1)) * moist(ngrass+fueltotal:,:,:,1)
    !    print*,'max value of moisture after adjustment = ',maxval(outarray%lh20(:,:,:)+moist(ngrass+fueltotal:,:,:,1))
    !  endif
    !  
!
    !end subroutine grassmoistureoverride  

end module support
  
  
  
  
  
  
  
  
  
  
  
  