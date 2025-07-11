!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! define_variables defines all variables, both constant and 
! user-defined, used throughout the program
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
! Grid variables
!-----------------------------------------------------------------------
subroutine define_grid_variables
  use grid_variables, only : nx,ny,nz,rhof,sizescale,moist,fueldepth, &
    zs,zheight,topofile,aa1,dx,dy,dz,nfuel
  use infile_variables, only : infuel,inx,iny,inz,idx,idy,idz,iintpr, &
    intopofile,izs,izheight,iaa1,ifuelin
  use fuel_variables, only : ntreefueltypes,istem,tfuelbins,ngrass, &
    ntspecies,ilitter
  use io_variables, only : workdir,filesep
  use constant_variables, only : tolerance
  implicit none
  
  ! Local Variables
  integer :: i,j,k
  integer :: ii,jj
  integer :: xbot,xtop,ybot,ytop
  real :: cells,xfrac,yfrac
  real,external :: zcart
  real,dimension(2) :: xcor,ycor
  
  ! Executable Code
  ntreefueltypes = istem*2+tfuelbins
  nfuel = infuel+ngrass+ntspecies*ntreefueltypes
  print*,ntspecies
  if(ilitter.gt.0) nfuel=nfuel+ntspecies*tfuelbins
  allocate(rhof(nfuel,nx,ny,nz)); rhof(:,:,:,:)=0.0
  allocate(sizescale(nfuel,nx,ny,nz)); sizescale(:,:,:,:)=0.0
  allocate(moist(nfuel,nx,ny,nz)); moist(:,:,:,:)=0.0
  allocate(fueldepth(nfuel,nx,ny,nz)); fueldepth(:,:,:,:)=0.0
  
  !---------------------------------------------------------------------
  ! Create topo layer (Should be adjusted for non-flat topo)
  !---------------------------------------------------------------------
  
  allocate(zs(nx,ny))
  allocate(zheight(nx,ny,nz))
  if (topofile.eq.'flat'.or.topofile.eq.'') then ! No topo
    zs(:,:)=0.0
    print *,'Not using target topo'
  else ! Normal topo
    print *,'Reading target topo file = ',topofile
    open (1,file=TRIM(TRIM(workdir)//filesep)//topofile, &
      form='unformatted',status='old')
    read (1) zs
    close (1)
  endif
  
  do i=1,nx
    do j=1,ny
      do k=1,nz
        if (aa1.eq.0) then
          zheight(i,j,k) = zs(i,j)+(k-1)*dz
        else
          zheight(i,j,k) = zcart(aa1,(k-1)*dz,nz,dz,zs(i,j))
        endif
        if (i.eq.1.and.j.eq.1)  &
          print*,'cell',k,'bottom height',zheight(i,j,k)
      enddo
    enddo
  enddo
  
  if(minval(zs).gt.0)then ! Reduce topo values to least common value
    zs  = zs-minval(zs)
    open (2,file=TRIM(TRIM(workdir)//filesep)//'toporeduced.dat', &
      form='unformatted',status='unknown')
    write(2) zs
    close(2)
  endif
  
  if(ifuelin.eq.1)then
    allocate(izs(inx,iny))
    allocate(izheight(inx,iny,inz))
  
    if(inx.ne.nx.or.iny.ne.ny.or.inz.ne.nz.or. &
      abs(idx-dx).lt.tolerance.or.abs(idy-dy).lt.tolerance.or. &
      abs(idz-dz).lt.tolerance.or.abs(iaa1-aa1).lt.tolerance.or. &
      topofile.ne.intopofile) &
      iintpr=1
  
    if (iintpr.eq.0) then
      izs(:,:)=zs(:,:)
      izheight(:,:,:)=zheight(:,:,:)
    else  ! Topo with existing fuels
      if (intopofile.eq.'flat'.or.intopofile.eq.'')then ! No previous topo
        izs(:,:)=0.0
        print *,'Not using previous topo'
      else  ! Normal previous topo
        print *,'Reading previous fuel topo file = ',intopofile
        open (2,file=TRIM(TRIM(workdir)//filesep)//intopofile, &
          form='unformatted',status='old')
        read (2) izs
        close (2)
      endif
      izs = izs-minval(zs)
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
    endif
  endif
end subroutine define_grid_variables

!-----------------------------------------------------------------------
! Variables unique to the fuels_create
!-----------------------------------------------------------------------
subroutine define_fuel_variables
  use grid_variables, only : nx,ny,nz
  use fuel_variables, only : igrass,grhof,ngrass,gsizescale,gmoist, &
    gfueldepth,ilitter,ntspecies,tfuelbins,lrhof,lsizescale,lmoist, &
    lfueldepth,trhof,tsizescale,tmoist,tfueldepth,itrees,trhofmax, &
    trhofmaxindex,t2bulkdensity,ntreefueltypes,ntrees
  use io_variables, only : verbose
  implicit none
  
  ! Local Variables
  integer :: ift,i,it
  real,allocatable :: averagedensity(:),trhofmaxtmp(:)
  
  ! Executable Code
  if (igrass.ne.0) then
    allocate(grhof(ngrass,nx,ny,nz))
    grhof(:,:,:,:)=0.0
    allocate(gsizescale(ngrass,nx,ny,nz)); gsizescale(:,:,:,:)=0.0
    allocate(gmoist(ngrass,nx,ny,nz)); gmoist(:,:,:,:)=0.0
    allocate(gfueldepth(ngrass,nx,ny)); gfueldepth(:,:,:)=0.0
  endif
  if (ilitter.ne.0) then
    it=ntspecies*tfuelbins
    allocate(lrhof(it,nx,ny,nz)); lrhof(:,:,:,:)=0.0
    allocate(lsizescale(it,nx,ny,nz)); lsizescale(:,:,:,:)=0.0
    allocate(lmoist(it,nx,ny,nz)); lmoist(:,:,:,:)=0.0
    allocate(lfueldepth(it,nx,ny)); lfueldepth(:,:,:)=0.0
  endif
  it=ntspecies*ntreefueltypes
  allocate(trhof(it,nx,ny,nz)); trhof(:,:,:,:)=0.0
  allocate(tsizescale(it,nx,ny,nz)); tsizescale(:,:,:,:)=0.0
  allocate(tmoist(it,nx,ny,nz)); tmoist(:,:,:,:)=0.0
  allocate(tfueldepth(it,nx,ny)); tfueldepth(:,:,:)=0.0
  
  !set name of output files
  if (verbose.eq.1) call define_newtreefile
  
  !---------------------------------------------------------------------
  ! Tree variables unique to the tree fuels_create
  !---------------------------------------------------------------------
  if(itrees.eq.2.or.itrees.eq.3) then
    call treelist_readin
  else if(itrees.eq.7) then
    call treelist_fastfuels
  endif
  
  ! Set maximum density tolerances for different species (unless values provided by user)
  if(itrees.gt.0)then
    allocate(trhofmaxindex(ntspecies))
    do ift=1,ntspecies
      if(trhofmax(ift).eq.0.) then
        allocate(averagedensity(ntrees(ift)))
        do i=1,ntrees(ift)
          averagedensity(i)=sum(t2bulkdensity(i,:,ift))
        enddo
        trhofmax(ift)=3./2.*maxval(averagedensity)
        deallocate(averagedensity)
      endif
    enddo
    allocate(trhofmaxtmp(ntspecies))
    trhofmaxtmp=trhofmax(:)
    do ift=1,ntspecies
      trhofmaxindex(ift)=maxloc(trhofmaxtmp,dim=1)
      trhofmaxtmp(maxloc(trhofmaxtmp,dim=1))=minval(trhofmax)-1
    enddo
    deallocate(trhofmaxtmp)
  endif
end subroutine define_fuel_variables

!-----------------------------------------------------------------------
!  
!-----------------------------------------------------------------------
subroutine define_newtreefile
  use fuel_variables, only: newtreefile, itrees, treefile
  implicit none  

  ! Local Variables
  integer(8) :: timeit
  character(len=30) :: dateit
  character(len=24) :: resultit
  integer :: n

  ! Executable Variables
  if (itrees.eq.0) then 
    newtreefile = '_blank'
  else 
    !------Open file for writing treelist to replicate files-------
    timeit = time8()
    call ctime(timeit,dateit)
    do n = 1, len(TRIM(dateit))
      resultit(n:n) = dateit(n:n)
      if ((dateit(n:n) == ' ').or.(dateit(n:n) == ':')) then
        resultit(n:n) = '_'
      endif
    enddo
    newtreefile = treefile(1:len_trim(treefile)-len_trim('.txt'))// &
      '_'//TRIM(resultit)
  endif
end subroutine define_newtreefile
