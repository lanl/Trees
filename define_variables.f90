!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! define_variables defines all variables, both constant and 
! user-defined, used throughout the program
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine define_constant_variables
!-----------------------------------------------------------------
! Constant variables
!-----------------------------------------------------------------
use constant_variables
implicit none

! Executable Code
PI = 3.14159

end subroutine define_constant_variables

subroutine define_grid_variables
!-----------------------------------------------------------------
! Grid variables
!-----------------------------------------------------------------
use grid_variables
use infile_variables
use baseline_variables

implicit none

! Local Variables
integer i,j,k
integer ii,jj
integer xbot,xtop,ybot,ytop
real    cells,xfrac,yfrac
real,external:: zcart
real,dimension(2):: xcor,ycor
     
! Executable Code
ntreefueltypes = istem*2+tfuelbins
nfuel = infuel+ngrass+ntspecies*ntreefueltypes
if(ilitter.gt.0) nfuel=nfuel+ntspecies*tfuelbins
allocate(rhof(nfuel,nx,ny,nz)); rhof(:,:,:,:)=0.0
allocate(sizescale(nfuel,nx,ny,nz)); sizescale(:,:,:,:)=0.0
allocate(moist(nfuel,nx,ny,nz)); moist(:,:,:,:)=0.0
allocate(fueldepth(nfuel,nx,ny,nz)); fueldepth(:,:,:,:)=0.0

!-----------------------------------------------------------------
! Create topo layer (Should be adjusted for non-flat topo)
!-----------------------------------------------------------------

allocate(zs(nx,ny))
allocate(izs(inx,iny))
allocate(zheight(nx,ny,nz))
allocate(izheight(inx,iny,inz))

if(ifuelin.eq.1.and.(inx.ne.nx.or.idx.ne.dx.or. &
  iny.ne.ny.or.idy.ne.dy.or.inz.ne.nz.or.idz.ne.dz &
  .or.aa1.ne.iaa1)) iintpr=1

if (topofile.eq.'flat'.or.topofile.eq.'') then ! No topo
  zs(:,:)=0.0
  izs(:,:)=0.0
  print *,'Not using topo'      
elseif (iintpr.eq.1) then ! Topo with existing fuels
  print *,'Reading given fuel topo file = ',topofile
  open (1,file=topofile,form='unformatted',status='old')
  read (1) izs
  close (1)
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
else ! Normal topo
  print *,'Reading topo file = ',topofile
  open (1,file=topofile,form='unformatted',status='old')
  read (1) zs
  close (1)
endif
if(minval(zs).gt.0)then ! Reduce topo values to least common value
  izs = izs-minval(zs)
  zs  = zs-minval(zs)
  open (2,file='toporeduced.dat',form='unformatted',status='unknown')
  write(2) zs
  close(2)
endif

do i=1,nx
  do j=1,ny
    do k=1,nz
      if (aa1.eq.0) then
        zheight(i,j,k) = zs(i,j)+(k-1)*dz
      else
        zheight(i,j,k) = zcart(aa1,(k-1)*dz,nz,dz,zs(i,j))
      endif
      if(i.eq.1.and.j.eq.1) print*,'cell',k,'bottom height',zheight(i,j,k)
    enddo
  enddo
enddo
if (iintpr.eq.1) then ! Topo with existing fuels
  print*,'Readin fuel grid heights'
  do i=1,inx
    do j=1,iny
      do k=1,inz
        if (iaa1.eq.0) then
          izheight(i,j,k) = izs(i,j)+(k-1)*idz
        else
          izheight(i,j,k) = zcart(iaa1,(k-1)*idz,inz,idz,izs(i,j))
        endif
        if(i.eq.1.and.j.eq.1) print*,'cell',k,'bottom height',izheight(i,j,k)
      enddo
    enddo
  enddo
endif

end subroutine define_grid_variables

subroutine define_baseline_variables
!-----------------------------------------------------------------
! Variables unique to the baseline
!-----------------------------------------------------------------
use grid_variables
use baseline_variables

implicit none

! Local Variables
integer ift

! Executable Code
if (igrass.ne.0) then
  allocate(grhof(ngrass,nx,ny,nz)); grhof(:,:,:,:)=0.0
  allocate(gsizescale(ngrass,nx,ny,nz)); gsizescale(:,:,:,:)=0.0
  allocate(gmoist(ngrass,nx,ny,nz)); gmoist(:,:,:,:)=0.0
  allocate(gfueldepth(ngrass,nx,ny)); gfueldepth(:,:,:)=0.0
endif
if (itrees.ne.0) then
  allocate(trhof(ntspecies*ntreefueltypes,nx,ny,nz)); trhof(:,:,:,:)=0.0
  allocate(tsizescale(ntspecies*ntreefueltypes,nx,ny,nz)); tsizescale(:,:,:,:)=0.0
  allocate(tmoist(ntspecies*ntreefueltypes,nx,ny,nz)); tmoist(:,:,:,:)=0.0
  allocate(tfueldepth(ntspecies*ntreefueltypes,nx,ny)); tfueldepth(:,:,:)=0.0
endif
if (ilitter.ne.0) then
  allocate(lrhof(ntspecies*tfuelbins,nx,ny,nz)); lrhof(:,:,:,:)=0.0
  allocate(lsizescale(ntspecies*tfuelbins,nx,ny,nz)); lsizescale(:,:,:,:)=0.0
  allocate(lmoist(ntspecies*tfuelbins,nx,ny,nz)); lmoist(:,:,:,:)=0.0
  allocate(lfueldepth(ntspecies*tfuelbins,nx,ny)); lfueldepth(:,:,:)=0.0
endif

!-----------------------------------------------------------------
! Groundfuel variables unique to the ground fuels baseline
!-----------------------------------------------------------------
if (igrass.eq.1) then
  allocate(grho(ngrass))
  allocate(gmoisture(ngrass))
  allocate(gss(ngrass))
  allocate(gdepth(ngrass))
  
  open (1,file=grassfile)
  read (1,*) grho       ! bulk density of grass [kg/m3]
  read (1,*) gmoisture  ! moisture content of grass
  read (1,*) gss        ! size scale of grass [m]
  read (1,*) gdepth     ! depth of grass [m]
  close (1)
elseif (igrass.eq.2) then
  print*,'Reading 2D ground fuel file'
  open (1,file=grassfile,form='unformatted',status='old')
  do ift=1,ngrass
    read (1) grhof(ift,:,:,1)      ! bulk density of grass [kg/m3]
    read (1) gmoist(ift,:,:,1)     ! moisture content of grass
    read (1) gsizescale(ift,:,:,1) ! size scale of grass [m]
    read (1) gfueldepth(ift,:,:) ! depth of grass [m]
  enddo
  close (1)
endif

!-----------------------------------------------------------------
! Tree variables unique to the tree baseline
!-----------------------------------------------------------------
if (itrees.eq.1) then
  call treesGeneral_readin
else if(itrees.eq.2.or.itrees.eq.3) then
  call treelist_readin
else if(itrees.eq.7) then
  call treelist_fastfuels
endif

if (ilitter.gt.0) then
  allocate(lrho(ntspecies*tfuelbins))
  allocate(lmoisture(ntspecies*tfuelbins))
  allocate(lss(ntspecies*tfuelbins))
  allocate(ldepth(ntspecies*tfuelbins))
  
  open (3,file=litterfile)
  read (3,*) lrho       ! bulk density of litter [kg/m3]
  read (3,*) lmoisture  ! moisture content of litter [fraction]
  read (3,*) lss        ! size scale of litter [m]
  read (3,*) ldepth     ! depth of litter [m]
  close (3)
endif

end subroutine define_baseline_variables

subroutine define_duet_variables
!-----------------------------------------------------------------
! Variables unique to the duet
!-----------------------------------------------------------------
use grid_variables, only : nx,ny,nz
use baseline_variables, only : ntspecies,tfuelbins
use duet_variables, only : windprofile,winddatafile,StepsPerYear, &
  YearsSinceBurn,uavg,vavg,VAR,Umean,Vmean,Uvar,Vvar,vterminal, &
  fuelSA,Froude,droptime,leafdropfreq,decay,speciesfile,lrhofT, &
  grhofT
implicit none

! Local Variables
integer :: yt,ift
integer :: fuelTotal,periodTotal
real :: fuelMass
real,dimension(6) :: temp_array

! Executable Code
periodTotal=YearsSinceBurn*StepsPerYear
if(windprofile.eq.0) then
  allocate(uavg(periodTotal),vavg(periodTotal),VAR(periodTotal,2))
  open (100,file=winddatafile)  ! wind information for area per year
  do yt=1,periodTotal
    read(100,*) temp_array
    uavg(yt)=temp_array(3)
    vavg(yt)=temp_array(4)
    VAR(yt,:)=temp_array(5:6)
  enddo
  close(100)
elseif(windprofile.eq.1)then
  allocate(Umean(nx,ny,nz),Vmean(nx,ny,nz),Uvar(nx,ny,nz),VVar(nx,ny,nz))
  open (105,file='Umean.dat', form='unformatted', status='old')
  read (105) Umean
  close(105)

  open (106,file='Vmean.dat', form='unformatted', status='old')
  read (106) Vmean
  close(106)

  open (107,file='Uvar.dat', form='unformatted', status='old')
  read (107) Uvar
  close(107)

  open (108,file='Vvar.dat', form='unformatted', status='old')
  read (108) Vvar
  close(108)
endif

! Species data
open (99,file=speciesfile)
fuelTotal=ntspecies*tfuelbins
allocate(fuelSA(fuelTotal),leafdropfreq(fuelTotal),decay(fuelTotal), &
  droptime(fuelTotal),vterminal(fuelTotal),Froude(fuelTotal))
do ift=1,ntspecies*tfuelbins
  read(99,*) fuelMass,fuelSA(ift),leafdropfreq(ift),decay(ift),droptime(ift)
  fuelMass=fuelMass/1000.
  fuelSA(ift)=fuelSA(ift)/10000.
  ! Terminal velocity
  vterminal(ift)=sqrt((2.*fuelMass*9.81)/(0.6*1.225*fuelSA(ift)))
  ! 9.81 is acceleration of gravity
  ! 0.6 is drag coefficient
  ! 1.225 is density of air at STP
  Froude(ift) = 0.01/(sqrt(9.81*sqrt(fuelSA(ift))))
enddo
close(99)

allocate(lrhofT(fuelTotal,nx,ny,periodTotal)); lrhofT=0.
allocate(grhofT(fuelTotal,nx,ny,periodTotal)); grhofT=0.

end subroutine define_duet_variables
