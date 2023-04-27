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
use baseline_variables
use duet_variables, only : inputprogram

implicit none
     
integer :: i,j,k


! Executable Code 

  ntreefueltypes = istem*2+tfuelbins
  !fueltotal = fueltotal !infuel+ntspecies*tfuelbins
  allocate(rhof(2*fueltotal+ngrass,nx,ny,nz)); rhof(:,:,:,:)=0.0
  allocate(sizescale(2*fueltotal+ngrass,nx,ny,nz)); sizescale(:,:,:,:)=0.0
  allocate(moist(2*fueltotal+ngrass,nx,ny,nz)); moist(:,:,:,:)=0.0
  allocate(zheight(nx,ny,nz)); zheight(:,:,:)=0.0

  allocate(trhof(fueltotal,nx,ny,nz)); trhof(:,:,:,:)=0.0
  allocate(tsizescale(fueltotal,nx,ny,nz)); tsizescale(:,:,:,:)=0.0
  allocate(tmoist(fueltotal,nx,ny,nz)); tmoist(:,:,:,:)=0.0
  allocate(tfueldepth(fueltotal,nx,ny,nz)); tfueldepth(:,:,:,:)=0.0


if(inputprogram.eq.1) then
  open(unit=100,file='fuelfiles_RMSZ.dat',form='unformatted',status='unknown')
    read(100) trhof(1:fueltotal,:,:,:),tmoist(1:fueltotal,:,:,:),&
    tsizescale(1:fueltotal,:,:,:),zheight(:,:,:)
  close(100)
  print*,'Input Program is Trees'
  print*,'minval, maxval of rhof = ',minval(trhof),maxval(trhof)
  print*,'minval, maxval of moist = ',minval(tmoist),maxval(tmoist)
  print*,'minval, maxval of sizescale = ',minval(tsizescale),maxval(tsizescale)
  print*,'minval, maxval of zheight = ',minval(zheight),maxval(zheight)
elseif(inputprogram.eq.2) then
  print*,'Input Program is FastFuels 3D'  
  do k=1,nz
    do j=1,ny
      do i=1,nx
        zheight(i,j,k) = (k-1)*dz
      enddo
    enddo
  enddo
  print*,'Max and Min of zheight = ',maxval(zheight),minval(zheight)
  if(any(isNaN(zheight))) print*,'zheight is NaN'
endif

!print*,'minval, maxval of rhof = ',minval(rhof),maxval(rhof)
!print*,'minval, maxval of moist = ',minval(moist),maxval(moist)
!print*,'minval, maxval of sizescale = ',minval(sizescale),maxval(sizescale)
!print*,'minval, maxval of zheight = ',minval(zheight),maxval(zheight)


end subroutine define_grid_variables

subroutine define_baseline_variables
!-----------------------------------------------------------------
! Variables unique to the baseline
!-----------------------------------------------------------------
use grid_variables, only : nx,ny,nz
use baseline_variables, only : fueltotal,lrhof,lsizescale,lmoist,lfueldepth

implicit none
! Executable Code

allocate(lrhof(fueltotal,nx,ny,nz)); lrhof(:,:,:,:)=0.0
allocate(lsizescale(fueltotal,nx,ny,nz)); lsizescale(:,:,:,:)=0.0
allocate(lmoist(fueltotal,nx,ny,nz)); lmoist(:,:,:,:)=0.0
allocate(lfueldepth(fueltotal,nx,ny)); lfueldepth(:,:,:)=0.0

end subroutine define_baseline_variables

subroutine define_species_variables
!-----------------------------------------------------------------
! Imports the species database for FIA codes
!-----------------------------------------------------------------

use species_variables

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

subroutine define_3Dspecies_variables

use species_variables

open(unit=5, file='FIA_FastFuels_fin_fulllist_populated.txt', &
  status='old', access='sequential', form='formatted')
  do i=1,290
    read(5,*) SPECINFO(i)
  enddo
close(5)

end subroutine define_3Dspecies_variables

subroutine define_duet_variables
!-----------------------------------------------------------------
! Variables unique to the duet
!-----------------------------------------------------------------
use grid_variables, only : nx,ny,nz
use baseline_variables, only : fueltotal,ngrass
use duet_variables, only : windprofile,StepsPerYear, &
  YearsSinceBurn,uavg,vavg,VAR,ustd,vstd,Umean,Vmean,Uvar,Vvar,vterminal, &
  fuelSA,Froude,droptime,leafdropfreq,decay,speciesfile,lrhofT, &
  grhofT,dragco,lafdT,gafdT,lmoistT,gmoistT,lssT,gssT,randomwinds, &
  moistspec,ssspec,compact,periodTotal
use species_variables
implicit none

! Local Variables
integer :: i,yt,ift,s,low,high,loc
real :: fuelMass,dragslope,dragb,a

! Executable Code
periodTotal=YearsSinceBurn*StepsPerYear
allocate(uavg(periodTotal),vavg(periodTotal))
allocate(ustd(periodTotal),vstd(periodTotal),VAR(periodTotal,2))
if(windprofile.eq.0) then
  print*,'!----!----!----!----!----!----!----!----!----!----!----!'
  print*,'Windprofile is taken from a user provided in fuellist'
  open(3,file='windprofile.dat',form='unformatted',status='old')
  read(3) uavg,vavg,ustd,vstd
  close(3)
  do yt=1,periodTotal
    VAR(yt,1)=ustd(yt)
    VAR(yt,2)=vstd(yt)
  enddo
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

allocate(fuelSA(fueltotal),leafdropfreq(fueltotal),decay(fueltotal), &
droptime(fueltotal),vterminal(fueltotal),Froude(fueltotal), &
moistspec(fueltotal),dragco(fueltotal),ssspec(fueltotal),compact(fueltotal))

!print*,'Works up to point 1'
! Species data
if (iFIA.eq.1) then
  allocate(FIA(fueltotal)); FIA(:) = 0
  !allocate(FIAtemp(fueltotal)); FIAtemp(:) = 0
  !print*,'Works up to point 2, FIA = ',FIA
  open (79,file='FIA.dat',form='unformatted',status='old',action='read',access='stream')
  !print*,'Works to open the file'
    do i=1,fueltotal
      read (79) FIA(i)
      !print*,'Works to read the file'
    enddo
  close(79)
  !print*,'Works up to point 3, FIAtemp = ',FIAtemp
  !FIA = int(FIAtemp)
  !print*,'Works up to point 3, FIA = ',FIA
  ift=1
  if(iFIAspecies.eq.0) then
    print*,'!----!----!----!----!----!----!----!----!----!----!----!----!----!----!'
    print*,'USING SPECIES GROUPS... check below to make sure you have inputted groups'
    print*,'FIA = ',FIA
    print*,'!----!----!----!----!----!----!----!----!----!----!----!----!----!----!'
    do ift=1,fueltotal
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
    do ift=1,fueltotal
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
  do ift=1,fueltotal
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

allocate(lrhofT(fueltotal,nx,ny,periodTotal)); lrhofT(:,:,:,:)=0.0
allocate(grhofT(ngrass,nx,ny,periodTotal)); grhofT(:,:,:,:)=0.0
allocate(lafdT(fueltotal,nx,ny,periodTotal)); lafdT(:,:,:,:)=0.0
allocate(gafdT(ngrass,nx,ny,periodTotal)); gafdT(:,:,:,:)=0.0
allocate(lmoistT(fueltotal,nx,ny,periodTotal)); lmoistT(:,:,:,:)=0.0
allocate(gmoistT(ngrass,nx,ny,periodTotal)); gmoistT(:,:,:,:)=0.0
allocate(lssT(fueltotal,nx,ny,periodTotal)); lssT(:,:,:,:)=0.0
allocate(gssT(ngrass,nx,ny,periodTotal)); gssT(:,:,:,:)=0.0

end subroutine define_duet_variables

subroutine define_3Dduet_variables
!-----------------------------------------------------------------
! Variables unique to the duet
!-----------------------------------------------------------------
use constant_variables
use grid_variables, only : nx,ny
use baseline_variables, only : fueltotal,ngrass
use duet_variables, only : StepsPerYear, &
  YearsSinceBurn,uavg,vavg,VAR,vterminal, &
  fuelSA,Froude,droptime,leafdropfreq,decay,lrhofT, &
  grhofT,dragco,lafdT,gafdT,lmoistT,gmoistT,lssT,gssT,randomwinds, &
  moistspec,ssspec,compact,periodTotal
use species_variables
use FF_variables, only : windvary,winddirection
implicit none

real :: a
integer :: low,high,yt

periodTotal=YearsSinceBurn*StepsPerYear
allocate(uavg(periodTotal),vavg(periodTotal))
allocate(VAR(periodTotal,2))

print*,'!----!----!----!----!----!----!----!----!----!----!----!' 
print*,'Windprofile is randomly generated averages over whole domain' 
print*,'multiplication factor for winds = ',randomwinds
low=1*randomwinds
high=3*randomwinds
do yt=1,periodTotal
  a = rand()
  a = 2.0*a - 1.0
  if(winddirection.ge.180) then
    a = abs(a)
  else
    a = -abs(a)
  endif

  if (a.gt.0) then
    uavg(yt) = ((a*real(low+high))-real(low))*cos(-(winddirection+90)*PI/180)
    vavg(yt) = ((a*real(low+high))-real(low))*sin(-(winddirection+90)*PI/180)
  else 
    uavg(yt) = ((a*real(low+high))+real(low))*cos(-(winddirection+90)*PI/180)
    vavg(yt) = ((a*real(low+high))+real(low))*sin(-(winddirection+90)*PI/180)
  endif
  if(windvary.eq.0) then
    a = rand()    
    a = (-(winddirection+windvary/2+90)*PI/180)*a
    if (a.gt.0) then
      VAR(yt,1) = 0.25*((a*real(low+high))-real(low))
    else
      VAR(yt,1) = 0.25*((a*real(low+high))+real(low))
    endif
    a = rand()     
    a = (-(winddirection-windvary/2+90)*PI/180)*a
    if (a.gt.0) then
      VAR(yt,2) = 0.25*((a*real(low+high))-real(low))
    else
      VAR(yt,2) = 0.25*((a*real(low+high))+real(low))
    endif
  else
    a = rand()    
    a = 2.0*a - 1.0
    if (a.gt.0) then
      VAR(yt,1) = 0.25*((a*real(low+high))-real(low))
    else
      VAR(yt,1) = 0.25*((a*real(low+high))+real(low))
    endif
    a = rand()     
    a = 2.0*a - 1.0
    if (a.gt.0) then
      VAR(yt,2) = 0.25*((a*real(low+high))-real(low))
    else
      VAR(yt,2) = 0.25*((a*real(low+high))+real(low))
    endif
  endif
enddo

print*,'uavg = ',uavg
print*,'vavg = ',vavg
print*,'uvar = ',VAR(:,1)
print*,'vvar = ',VAR(:,2)
print*,'Each column is a year'
print*,'!----!----!----!----!----!----!----!----!----!----!----!' 

allocate(fuelSA(fueltotal),leafdropfreq(fueltotal),decay(fueltotal), &
droptime(fueltotal),vterminal(fueltotal),Froude(fueltotal), &
moistspec(fueltotal),dragco(fueltotal),ssspec(fueltotal),compact(fueltotal))

allocate(lrhofT(fueltotal,nx,ny,periodTotal)); lrhofT(:,:,:,:)=0.0
allocate(grhofT(ngrass,nx,ny,periodTotal)); grhofT(:,:,:,:)=0.0
allocate(lafdT(fueltotal,nx,ny,periodTotal)); lafdT(:,:,:,:)=0.0
allocate(gafdT(ngrass,nx,ny,periodTotal)); gafdT(:,:,:,:)=0.0
allocate(lmoistT(fueltotal,nx,ny,periodTotal)); lmoistT(:,:,:,:)=0.0
allocate(gmoistT(ngrass,nx,ny,periodTotal)); gmoistT(:,:,:,:)=0.0
allocate(lssT(fueltotal,nx,ny,periodTotal)); lssT(:,:,:,:)=0.0
allocate(gssT(ngrass,nx,ny,periodTotal)); gssT(:,:,:,:)=0.0


end subroutine define_3Dduet_variables

subroutine define_ff_variables

use grid_variables
use FF_variables

allocate(FFrhof(nx,ny,nz),FFmoist(nx,ny,nz))
allocate(FFspec(nx,ny,nz))
allocate(surfrhof(nx,ny),surfdepth(nx,ny))

surfrhof = 0.0
surfdepth = 0.0

open(1,file='treesrhof.dat',form='unformatted',status='old')
read(1) FFrhof
close(1)

open(2,file='treesmoist.dat',form='unformatted',status='old')
read(2) FFmoist
close(2)

open(unit=1,file='treesspcd.dat',form='unformatted',status='unknown')
read(1) FFSpec
close(1)

!if(any(isNaN(FFSpec))) print*,'Reading in NaNs for FF Data species'

print*,'Shape of FFrhof = ',shape(FFrhof)
print*,'Max,Min of FFrhof = ',maxval(FFrhof),minval(FFrhof)
if(any(isNaN(FFrhof))) print*,'Reading in NaNs for FF Data rhof'
print*,'Shape of FFmoist = ',shape(FFmoist)
print*,'Max,Min of FFmoist = ',maxval(FFmoist),minval(FFmoist)
if(any(isNaN(FFmoist))) print*,'Reading in NaNs for FF Data moisture'
print*,'Shape of FFspec = ',shape(FFspec)
print*,'Max,Min of FFspec = ',maxval(FFspec),minval(FFspec)

end subroutine define_ff_variables
