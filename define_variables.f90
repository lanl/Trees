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
!-----------------------------------------------------------------
subroutine define_constant_variables
!-----------------------------------------------------------------
! Constant variables
!-----------------------------------------------------------------
use constant_variables
implicit none

! Executable Code
PI = 3.14159 !YUM

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
      .or.aa1.ne.iaa1.or.topofile.ne.intopofile)) &
      iintpr=1

if (iintpr.eq.0) then
  if (topofile.eq.'flat'.or.topofile.eq.'') then ! No topo
    zs(:,:)=0.0
    izs(:,:)=0.0
    print *,'Not using target topo'
  else ! Normal topo
    print *,'Reading target topo file = ',topofile
    open (1,file=topofile,form='unformatted',status='old')
    read (1) zs
    close (1)
    izs(:,:)=zs(:,:)
  endif
endif

if (iintpr.eq.1) then ! Topo with existing fuels
  if (topofile.eq.'flat'.or.topofile.eq.'')then ! No target topo JO
    zs(:,:)=0.0
    print *,'Not using target topo'
  else  ! Normal target topo
    print *,'Define Varibles: Reading target topo file = ',topofile
    open (1,file=topofile,form='unformatted',status='old')
    read (1) zs(:,:)
    close (1)
    !izs(:,:)=zs(:,:)
  endif
  if (intopofile.eq.'flat'.or.intopofile.eq.'')then ! No previous topo
    izs(:,:)=0.0
    print *,'Not using previous topo'
  else  ! Normal previous topo
    print *,'Reading previous fuel topo file = ',intopofile    !JO
    open (2,file=intopofile,form='unformatted',status='old') !JO
    read (2) izs
    close (2)
  endif
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
else
  izheight(:,:,:)=zheight(:,:,:)
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
  allocate(grhof(ngrass,nx,ny,nz))
  grhof(:,:,:,:)=0.0
  allocate(gsizescale(ngrass,nx,ny,nz)); gsizescale(:,:,:,:)=0.0
  allocate(gmoist(ngrass,nx,ny,nz)); gmoist(:,:,:,:)=0.0
  allocate(gfueldepth(ngrass,nx,ny)); gfueldepth(:,:,:)=0.0
endif
if (ilitter.ne.0) then
  allocate(lrhof(ntspecies*tfuelbins,nx,ny,nz)); lrhof(:,:,:,:)=0.0
  allocate(lsizescale(ntspecies*tfuelbins,nx,ny,nz)); lsizescale(:,:,:,:)=0.0
  allocate(lmoist(ntspecies*tfuelbins,nx,ny,nz)); lmoist(:,:,:,:)=0.0
  allocate(lfueldepth(ntspecies*tfuelbins,nx,ny)); lfueldepth(:,:,:)=0.0
endif
allocate(trhof(ntspecies*ntreefueltypes,nx,ny,nz)); trhof(:,:,:,:)=0.0
allocate(tsizescale(ntspecies*ntreefueltypes,nx,ny,nz)); tsizescale(:,:,:,:)=0.0
allocate(tmoist(ntspecies*ntreefueltypes,nx,ny,nz)); tmoist(:,:,:,:)=0.0
allocate(tfueldepth(ntspecies*ntreefueltypes,nx,ny)); tfueldepth(:,:,:)=0.0


!-----------------------------------------------------------------
! Groundfuel variables unique to the ground fuels baseline
!-----------------------------------------------------------------

! Grass and litter inputs defined in io.f90 from fuellist

if (igrass.eq.2) then
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
end subroutine define_baseline_variables

subroutine define_species_variables
!-----------------------------------------------------------------
! Imports the species database for FIA codes
!-----------------------------------------------------------------

use species_variables
use baseline_variables, only:ntspecies


integer :: i,j,ct,fdrop=0,fstep=0
real :: fsa=0,fdecay=0,fmoist=0,fss=0
real :: fdrag=0,fvterm=0,ffrou=0,fcomp=0
open(unit=5, file='FIA_FastFuels_fin_fulllist_populated.txt', &
  status='old', access='sequential', form='formatted')
  do i=1,290
    read(5,*) SPECINFO(i)
  enddo
close(5)
if (iFIAspecies.eq.0) then
  do j=1,maxval(SPECINFO%sp_grp_num)
    ct=0
    fsa=0
    fdrop=0
    fdecay=0
    fmoist=0
    fstep=0
    fdrag=0
    fvterm=0
    ffrou=0
    fcomp=0
    fss=0
    do i=1,290
      if(SPECINFO(i)%sp_grp_num.eq.j) then
        ct=ct+1
        fsa=fsa+SPECINFO(i)%surfarea
        fdrop=fdrop+SPECINFO(i)%dropperyear
        fdecay=fdecay+SPECINFO(i)%decay
        fmoist=fmoist+SPECINFO(i)%moist
        fstep=fstep+SPECINFO(i)%stepperyear
        fdrag=fdrag+SPECINFO(i)%dragco
        fvterm=fvterm+SPECINFO(i)%vterminal
        ffrou=ffrou+SPECINFO(i)%froude
        fcomp=fcomp+SPECINFO(i)%compact
        fss=fss+SPECINFO(i)%sizescale
      endif
    enddo
    SPECgroups(j)%surfarea=fsa/real(ct)
    SPECgroups(j)%dropperyear=int(fdrop/ct)
    SPECgroups(j)%decay=fdecay/real(ct)
    SPECgroups(j)%moist=fmoist/real(ct)
    SPECgroups(j)%stepperyear=int(fstep/ct)
    SPECgroups(j)%dragco=fdrag/real(ct)
    SPECgroups(j)%vterminal=fvterm/real(ct)
    SPECgroups(j)%froude=ffrou/real(ct)
    SPECgroups(j)%compact=fcomp/real(ct)
    SPECgroups(j)%sizescale=fss/real(ct)
  enddo
endif

end subroutine define_species_variables

subroutine define_duet_variables
!-----------------------------------------------------------------
! Variables unique to the duet
!-----------------------------------------------------------------
use grid_variables, only : nx,ny,nz
use baseline_variables, only : ntspecies,tfuelbins
use duet_variables, only : windprofile,winddatafile,StepsPerYear, &
  YearsSinceBurn,uavg,vavg,VAR,ustd,vstd,Umean,Vmean,Uvar,Vvar,vterminal, &
  fuelSA,Froude,droptime,leafdropfreq,decay,speciesfile,lrhofT, &
  grhofT,dragco,lafdT,gafdT,lmoistT,gmoistT,lssT,gssT,randomwinds, &
  moistspec,ssspec,compact,periodTotal
use species_variables
implicit none

! Local Variables
integer :: yt,ift,s,low,high,loc
integer :: fuelTotal
real :: fuelMass,dragslope,dragb,a,monstep
real,dimension(6) :: temp_array 

! Executable Code
periodTotal=YearsSinceBurn*StepsPerYear
if(windprofile.eq.0) then
  print*,'!----!----!----!----!----!----!----!----!----!----!----!'
  print*,'Windprofile is taken from a user provided in fuellist'
  allocate(VAR(periodTotal,2))
  !open (100,file=winddatafile)  ! wind information for area per year
  do yt=1,periodTotal
    !read(100,*) temp_array
    !uavg(yt)=temp_array(3)
    !vavg(yt)=temp_array(4)
    VAR(yt,1)=ustd(yt)
    VAR(yt,2)=vstd(yt)
  enddo
  close(100)
elseif(windprofile.eq.1)then
  print*,'!----!----!----!----!----!----!----!----!----!----!----!'
  print*,'Windprofile is generated from higrad and averaged within each cell'
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

elseif(windprofile.eq.2) then
  print*,'!----!----!----!----!----!----!----!----!----!----!----!' 
  print*,'Windprofile is randomly generated averages over whole domain' 
  print*,'multiplication factor for winds = ',randomwinds
  low=1*randomwinds
  high=3*randomwinds

  allocate(uavg(periodTotal),vavg(periodTotal),VAR(periodTotal,2))

  do yt=1,periodTotal
    a = rand()
    a = 2.0*a - 1.0
    if (a.gt.0) then
      uavg(yt) = (a*real(low+high))-real(low)
    else 
      uavg(yt) = (a*real(low+high))+real(low)
    endif
    a = rand()
    a = 2.0*a - 1.0
    if (a.gt.0) then
      vavg(yt) = (a*real(low+high))-real(low)
    else
      vavg(yt) = (a*real(low+high))+real(low)
    endif
    a = rand()    
    a = 2.0*a - 1.0
    if (a.gt.0) then
      VAR(yt,1) = (a*real(low+high))-real(low)
    else
      VAR(yt,1) = (a*real(low+high))+real(low)
    endif
    a = rand()     
    a = 2.0*a - 1.0
    if (a.gt.0) then
      VAR(yt,2) = (a*real(low+high))-real(low)
    else
      VAR(yt,2) = (a*real(low+high))+real(low)
    endif
  enddo

endif
print*,'uavg = ',uavg
print*,'vavg = ',vavg
print*,'uvar = ',VAR(:,1)
print*,'vvar = ',VAR(:,2)
print*,'Each column is a year'
print*,'!----!----!----!----!----!----!----!----!----!----!----!' 

fuelTotal=ntspecies*tfuelbins
allocate(fuelSA(fuelTotal),leafdropfreq(fuelTotal),decay(fuelTotal), &
droptime(fuelTotal),vterminal(fuelTotal),Froude(fuelTotal), &
moistspec(fuelTotal),dragco(fuelTotal),ssspec(fuelTotal),compact(fuelTotal))

! Species data
if (iFIA.eq.1) then
  ift=1
  if(iFIAspecies.eq.0) then
    print*,'!----!----!----!----!----!----!----!----!----!----!----!----!----!----!'
    print*,'USING SPECIES GROUPS... check below to make sure you have inputted groups'
    print*,'FIA = ',FIA
    print*,'!----!----!----!----!----!----!----!----!----!----!----!----!----!----!'
    do ift=1,ntspecies*tfuelbins
      do s=1,10
        if(s.eq.FIA(ift)) then
          fuelSA(ift)=SPECgroups(s)%surfarea
          leafdropfreq(ift)=SPECgroups(s)%dropperyear
          decay(ift)=SPECgroups(s)%decay
          droptime(ift)=SPECgroups(s)%stepperyear
          vterminal(ift)=SPECgroups(s)%vterminal
          Froude(ift)=SPECgroups(s)%froude
          dragco(ift)=SPECgroups(s)%dragco
          moistspec(ift)=SPECgroups(s)%moist
          ssspec(ift)=SPECgroups(s)%sizescale
          compact(ift)=SPECgroups(s)%compact
          loc=s
          exit
        endif
      enddo
    enddo
  elseif (iFIAspecies.eq.1) then
    print*,'!----!----!----!----!----!----!----!----!----!----!----!----!----!----!'
    print*,'USING SPECIFIC SPECIES CODES... check below to make sure you have inputted species'
    print*,'FIA = ',FIA
    print*,'!----!----!----!----!----!----!----!----!----!----!----!----!----!----!'
    do ift=1,ntspecies*tfuelbins
      do s=1,290
        loc=0
        if(SPECINFO(s)%FIA_code.eq.FIA(ift)) then
          fuelSA(ift)=SPECINFO(s)%surfarea
          leafdropfreq(ift)=SPECINFO(s)%dropperyear
          decay(ift)=SPECINFO(s)%decay
          droptime(ift)=SPECINFO(s)%stepperyear
          vterminal(ift)=SPECINFO(s)%vterminal
          Froude(ift)=SPECINFO(s)%froude
          dragco(ift)=SPECINFO(s)%dragco
          moistspec(ift)=SPECINFO(s)%moist
          ssspec(ift)=SPECINFO(s)%sizescale
          compact(ift)=SPECINFO(s)%compact
          loc=s
          exit
        endif
      enddo
      !print*,'ift = ',ift
      !print*,'loc = ',loc
      !print*,'fuelSA(ift) = ',fuelSA(ift-1)
      !print*,'SPECINFO(s)%surfarea = ',SPECINFO(loc)%surfarea
    enddo
  endif  
elseif (iFIA.eq.0) then
  print*,'!----!----!----!----!----!----!----!----!----!----!----!----!----!----!'
  print*,'YOU ARE NOT USING THE FIA DATABASE - SPECIES FILE PROVIDED SEPARATELY'
  print*,'!----!----!----!----!----!----!----!----!----!----!----!----!----!----!'

    dragslope = 1.4/(sqrt(70.0) - 1.0)
    dragb = 0.6 - dragslope
  open (99,file=speciesfile)
  do ift=1,ntspecies*tfuelbins
    read(99,*) fuelMass,fuelSA(ift),leafdropfreq(ift),decay(ift),droptime(ift)
    fuelMass=fuelMass/1000.
    fuelSA(ift)=fuelSA(ift)/10000.
    dragco(ift)=dragslope*fuelSA(ift) + dragb !JSM
    ! Terminal velocity
    vterminal(ift)=sqrt((2.*fuelMass*9.81)/(dragco(ift)*1.225*fuelSA(ift)))
    ! 9.81 is acceleration of gravity
    ! 1.225 is density of air at STP
    Froude(ift) = 0.01/(sqrt(9.81*sqrt(fuelSA(ift))))
    moistspec(ift)=100
    compact(ift) = 0.5
  enddo
  close(99)
endif

droptime = ceiling(droptime/(12/StepsPerYear))

!print*,'Variable definitions:'
!print*,'fuelSA = ',fuelSA
!print*,'leafdropfreq = ',leafdropfreq
!print*,'decay = ',decay
!print*,'droptime = ',droptime
!print*,'vterminal = ',vterminal
!print*,'Froude = ',Froude
!print*,'dragco = ',dragco
!print*,'moistspec = ',moistspec
!print*,'ssspec = ',ssspec

allocate(lrhofT(fuelTotal,nx,ny,periodTotal)); lrhofT(:,:,:,:)=0.0
allocate(grhofT(fuelTotal,nx,ny,periodTotal)); grhofT(:,:,:,:)=0.0
allocate(lafdT(fuelTotal,nx,ny,periodTotal)); lafdT(:,:,:,:)=0.0
allocate(gafdT(fuelTotal,nx,ny,periodTotal)); gafdT(:,:,:,:)=0.0
allocate(lmoistT(fuelTotal,nx,ny,periodTotal)); lmoistT(:,:,:,:)=0.0
allocate(gmoistT(fuelTotal,nx,ny,periodTotal)); gmoistT(:,:,:,:)=0.0
allocate(lssT(fuelTotal,nx,ny,periodTotal)); lssT(:,:,:,:)=0.0
allocate(gssT(fuelTotal,nx,ny,periodTotal)); gssT(:,:,:,:)=0.0

end subroutine define_duet_variables
