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
    zs,zheight,topofile,aa1,dz,nfuel
  use infile_variables, only : infuel
  use fuel_variables, only : ntreefueltypes,istem,tfuelbins,ngrass, &
    ntspecies,ilitter
  use io_variables, only : workdir,filesep
  use constant_variables, only : tolerance
  implicit none
  
  ! Local Variables
  integer :: i,j,k
  real,external :: zcart
  
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
  allocate(zheight(nx,ny,nz+1))
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
      do k=1,nz+1
        if (aa1.lt.tolerance) then
          zheight(i,j,k) = zs(i,j)+(k-1)*dz
        else
          zheight(i,j,k) = zcart(aa1,(k-1)*dz,nz,dz,zs(i,j))
        endif
        if (i.eq.1.and.j.eq.1)  &
          print*,'cell',k,'bottom height',zheight(i,j,k)
      enddo
    enddo
  enddo
  
  if(abs(minval(zs)).lt.tolerance)then ! Reduce topo values to least common value
    zs  = zs-minval(zs)
    open (2,file=TRIM(TRIM(workdir)//filesep)//'toporeduced.dat', &
      form='unformatted',status='unknown')
    write(2) zs
    close(2)
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
  use constant_variables, only : tolerance
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
  
  !---------------------------------------------------------------------
  ! Tree variables unique to the tree fuels_create
  !---------------------------------------------------------------------
  if(itrees.gt.0) call treelist_readin
  
  ! Set maximum density tolerances for different species (unless values provided by user)
  if(itrees.gt.0)then
    allocate(trhofmaxindex(ntspecies))
    do ift=1,ntspecies
      if(trhofmax(ift).lt.tolerance) then
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
! define_newtreefile creates and new treelist with the verbose option
!-----------------------------------------------------------------------
subroutine define_newtreefile
  use fuel_variables, only: newtreefile,treefile
  implicit none  

  ! Local Variables
  character(len=30) :: dateit
  character(len=24) :: resultit
  integer :: n

  ! Executable Variables
  !------Open file for writing treelist to replicate files-------
  call ctime(time8(),dateit)
  do n = 1, len(TRIM(dateit))
    resultit(n:n) = dateit(n:n)
    if ((dateit(n:n) == ' ').or.(dateit(n:n) == ':')) then
      resultit(n:n) = '_'
    endif
  enddo
  newtreefile = treefile(1:len_trim(treefile)-len_trim('.txt'))// &
    '_'//TRIM(resultit)
end subroutine define_newtreefile
