!-----------------------------------------------------------------------
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
! grid_readin is a function which reads in .dat files for use in 
! FIRETEC or QUIC-Fire
!-----------------------------------------------------------------------
subroutine grid_readin
  use grid_variables, only : nx,ny,nz,dx,dy,dz,zs,zheight,rhof,moist, &
    sizescale,fueldepth,aa1,topofile
  use infile_variables, only : rhoffile,moistfile,ssfile,afdfile, &
    irhof,imoist,iss,iafd,infuel,inx,iny,inz,iintpr,idx,idy,idz, &
    izs,izheight,iaa1,intopofile
  use io_variables, only : filesep,workdir
  use constant_variables, only : tolerance
  implicit none
  
  ! Local Variables
  integer :: ift,i,j,k
  integer :: ii,jj,kk
  integer :: xbot,xtop,ybot,ytop,zbot,ztop
  real :: cells,xfrac,yfrac,zfrac
  real :: rhoftemp 
  real :: target_mass,actual_mass
  real,dimension(2) :: xcor,ycor
  real,external :: zcart
  
  ! Executable Code
  allocate(irhof(infuel,inx,iny,inz))
  allocate(imoist(infuel,inx,iny,inz))
  allocate(iss(infuel,inx,iny,inz))
  allocate(iafd(infuel,inx,iny,inz))
  allocate(izs(inx,iny))
  allocate(izheight(inx,iny,inz))
  
  print*,"Reading existing fuel files"
  ! Read in rhof
  open(unit=1,file=TRIM(TRIM(workdir)//filesep)//rhoffile, &
    form='unformatted',status='unknown')
  do ift=1,infuel
    read(1) irhof(ift,:,:,:)
  enddo
  close(1)
  irhof=max(0.,irhof)
  print*,'irhof: min=',minval(irhof),'max=',maxval(irhof)

  ! Read in moist
  open(unit=2,file=TRIM(TRIM(workdir)//filesep)//moistfile, &
    form='unformatted',status='unknown')
  do ift=1,infuel
    read(2) imoist(ift,:,:,:)
  enddo
  close(2)
  imoist=max(0.,imoist)
  print*,'imoist: min=',minval(imoist),'max=',maxval(imoist)

  ! Read in ss
  open(unit=3,file=TRIM(TRIM(workdir)//filesep)//ssfile, &
    form='unformatted',status='unknown')
  do ift=1,infuel
    read(3) iss(ift,:,:,:)
  enddo
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
  iss=max(0.,iss)
  print*,'iss: min=',minval(iss),'max=',maxval(iss)

  ! Read in afd
  open(unit=4,file=TRIM(TRIM(workdir)//filesep)//afdfile, &
    form='unformatted',status='unknown')
  do ift=1,infuel
    read(4) iafd(ift,:,:,:)
  enddo
  close(4)
  iafd=max(0.,iafd)
  print*,'iafd: min=',minval(iafd),'max=',maxval(iafd)
  
  ! Set grid characteristics  
  print*,'Readin topography'
  if (intopofile.eq.'flat'.or.intopofile.eq.'')then ! No previous topo
    izs(:,:)=0.0
  else  ! Normal previous topo
    print *,'Reading previous fuel topo file = ',intopofile
    open (2,file=TRIM(TRIM(workdir)//filesep)//intopofile, &
      form='unformatted',status='old')
    read (2) izs
    close (2)
  endif
  izs = izs-minval(izs)

  print*,'Readin fuel grid heights'
  do i=1,inx
    do j=1,iny
      do k=1,inz
        if (iaa1.eq.0) then
          izheight(i,j,k) = izs(i,j)+(k-1)*idz
        else
          izheight(i,j,k) = zcart(iaa1,(k-1)*idz,inz,idz,izs(i,j))
        endif
        if(i.eq.1.and.j.eq.1) print*,'cell',k,'bottom height', &
          izheight(i,j,k)
      enddo
    enddo
  enddo
  
  ! Interpolate read file onto FIRETEC grid
  if(inx.ne.nx.or.iny.ne.ny.or.inz.ne.nz.or. &
    abs(idx-dx).lt.tolerance.or.abs(idy-dy).lt.tolerance.or. &
    abs(idz-dz).lt.tolerance.or.abs(iaa1-aa1).lt.tolerance.or. &
    topofile.ne.intopofile) &
    iintpr=1
  
  if (iintpr.eq.1) then
    print*,"Interpolating readin fuel files to desired FIRETEC grid"
    do i=1,nx
      xcor(1) = (i-1)*dx ! Real x lower edge
      xcor(2) = i*dx     ! Real x upper edge
      xbot    = min(inx,floor(xcor(1)/idx+1)) ! Fuel readin grid x lower edge
      xtop    = min(inx,floor(xcor(2)/idx+1)) ! Fuel readin grid x upper edge
      do j=1,ny
        ycor(1) = (j-1)*dy ! Real y lower edge
        ycor(2) = j*dy     ! Real y upper edge
        ybot    = min(iny,floor(ycor(1)/idy+1)) ! Fuel readin grid y lower edge
        ytop    = min(iny,floor(ycor(2)/idy+1)) ! Fuel readin grid y upper edge 
        do k=1,nz-1
          zbot=inz
          ztop=inz
          do kk=1,inz
            if((izheight(ii,jj,kk+1)-izs(ii,jj)).gt. &
              (zheight(i,j,k)-zs(i,j)))then
              zbot=kk
              exit
            endif
          enddo
          do kk=zbot,inz
            if((izheight(ii,jj,kk)-izs(ii,jj)).ge. &
              (zheight(i,j,k+1)-zs(i,j)))then
              ztop=kk
              exit
            endif
          enddo
          cells = 0.
          do ii=xbot,xtop
            xfrac = (min(ii*idx,xcor(2))-max((ii-1)*idx,xcor(1)))/idx
            do jj=ybot,ytop
              yfrac = (min(jj*idy,ycor(2))-max((jj-1)*idy,ycor(1)))/idy
              do kk=zbot,ztop
                zfrac = min((izheight(ii,jj,kk+1)-izs(ii,jj)), &
                  (zheight(i,j,k+1)-zs(i,j))) - &
                  max((izheight(ii,jj,kk)-izs(ii,jj)) ,  &
                  (zheight(i,j,k)-zs(i,j)))/ &
                  (izheight(ii,jj,kk+1)-izheight(ii,jj,kk))
                cells = cells+xfrac*yfrac*zfrac
                do ift=1,infuel
                  rhoftemp=irhof(ift,ii,jj,kk)*xfrac*yfrac*zfrac
                  if(rhoftemp.gt.0)then
                    sizescale(ift,i,j,k) =  &
                      (rhof(ift,i,j,k)*sizescale(ift,i,j,k)+ &
                       rhoftemp*iss(ift,ii,jj,kk))/ &
                      (rhof(ift,i,j,k)+rhoftemp)
                    moist(ift,i,j,k) =  &
                      (rhof(ift,i,j,k)*moist(ift,i,j,k)+ &
                       rhoftemp*imoist(ift,ii,jj,kk))/ &
                      (rhof(ift,i,j,k)+rhoftemp)
                    fueldepth(ift,i,j,k) = &
                      (rhof(ift,i,j,k)*fueldepth(ift,i,j,k)+ &
                       rhoftemp*iafd(ift,ii,jj,kk))/ &
                      (rhof(ift,i,j,k)+rhoftemp)
                    rhof(ift,i,j,k)=rhof(ift,i,j,k)+rhoftemp
                  endif
                enddo
              enddo
            enddo
          enddo
          do ift=1,infuel
            if(cells.gt.0) rhof(ift,i,j,k)=rhof(ift,i,j,k)/cells
          enddo
        enddo
      enddo
    enddo
  else
    do ift=1,infuel
      rhof(ift,:,:,:)=irhof(ift,:,:,:)
      sizescale(ift,:,:,:)=iss(ift,:,:,:)
      moist(ift,:,:,:)=imoist(ift,:,:,:)
      fueldepth(ift,:,:,:)=iafd(ift,:,:,:)
    enddo
  endif
  
  ! Print out the target and actual fuel masses for comparisons sake
  target_mass = 0
  do ift=1,infuel
    do i=1,min(inx,nint(nx/dx*idx))
      do j=1,min(iny,nint(ny/dy*idy))
        do k=1,inz-1
          target_mass=target_mass+ &
            irhof(ift,i,j,k)*idx*idy*(izheight(i,j,k+1)-izheight(i,j,k))
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
          actual_mass=actual_mass+ &
            rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
        enddo
      enddo
    enddo
  enddo
  print*,'Readin actual fuel mass:',actual_mass
  print*,'Readin error:',actual_mass/target_mass*100,'%'
   
  deallocate(irhof,imoist,iss,iafd,izs,izheight)
end subroutine grid_readin
