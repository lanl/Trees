!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine will take a given tree canopy field and produce 
! surface fuels. It is an autonomous program that calculates how the 
! litter will fall from each tree, gather on the ground, and decay and 
! compact over time. Then it calculates where the grass can grow based 
! on the density of litter in that area and how much sunlight the area 
! will get.
! 
! Author: Jenna McDanold (2/21)
! Last Modified: 4/22 (AJJ)
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
subroutine Duet_Inputs
!-----------------------------------------------------------------
! Using canopy small cell info and grid location, calculate fall 
! time, then x and y displacement based on average wind info. Then
! using this info, find the location on the ground for the 
! displacement and horizontal and vertical stretch to place an 
! elliptical area for dispersing the fuel. Makes an elliptical 
! cylinder with volume 1, and find the area of that ellipse within 
! each grod cell containing the ellipse. This becomes the 
! percentage of the litter within the canopy cell that will be 
! placed within that higrad cell. 
!-----------------------------------------------------------------
use constant_variables, only : PI
use grid_variables, only : nx,ny,nz,dx,dy,zheight,zmax,rhof, &
  sizescale,moist,fueldepth
use infile_variables, only : infuel,ifuelin
use io_variables, only : controlseed,seedchange,singlefuel
use baseline_variables, only : ntspecies,tfuelbins,trhof, &
  tsizescale,tmoist,tfueldepth,grassconstant,litterconstant, &
  ngrass,gdepth,gmoisture,gss,grho,gmoistoverride,itrees,duet_ngrass
use duet_variables, only : vterminal,StepsPerYear,YearsSinceBurn, &
  Froude,droptime,windprofile,lrhofT,leafdropfreq,decay,grhofT, &
  uavg,vavg,VAR,ustd,vstd,Umean,Vmean,Uvar,Vvar,fuelSA,lmoistT,gmoistT, &
  lssT,gssT,lafdT,gafdT,compact,moistspec,ssspec,dragco,relhum, &
  periodTotal,litout,grassstep,randomwinds,inputprogram
use species_variables

implicit none

! Local variables
real,allocatable :: Drhof(:,:,:,:),Dss(:,:,:,:),Dmoist(:,:,:,:)
integer,allocatable ::a(:)
integer :: fueltotal
character(len=100) :: FuelFile,grassfile

! Executable code
!call define_duet_variables

fueltotal = infuel+ntspecies*tfuelbins
!print*,'INTS WE CARE ABOUT'
!print*,infuel,ntspecies,tfuelbins
!print*,fueltotal

allocate(Drhof(fueltotal,nx,ny,nz)); Drhof=0.0
allocate(Dss(fueltotal,nx,ny,nz)); DSS=0.0
allocate(Dmoist(fueltotal,nx,ny,nz)); Dmoist=0.0
allocate(a(4));a=1

!print*,'ARRAYS BEFORE FILL'
!print*,shape(Drhof)
!print*,shape(Dss)
!print*,shape(Dmoist)

if(ifuelin.eq.0.and.itrees.ne.1) then
  Drhof  = trhof
  Dss    = tsizescale
  Dmoist = tmoist
elseif(ifuelin.eq.0.and.itrees.eq.1) then
  Drhof  = rhof(1+ngrass:ngrass+fueltotal,:,:,:)
  Dss    = sizescale(1+ngrass:ngrass+fueltotal,:,:,:)
  Dmoist = moist(1+ngrass:ngrass+fueltotal,:,:,:)
elseif(ifuelin.eq.1.and.itrees.eq.0) then
  Drhof  = rhof(1:fueltotal,:,:,:)
  Dss    = sizescale(1:fueltotal,:,:,:)
  Dmoist = moist(1:fueltotal,:,:,:)  
else
  Drhof  = rhof
  Dss    = sizescale
  Dmoist = moist
  a = shape(rhof)
  fueltotal = a(1)
endif

!print*,'ARRAYS AFTER FILL'
!print*,Nx,Ny,Nz
!print*,shape(Drhof)
!print*,shape(Dss)
!print*,shape(Dmoist)
!print*,shape(rhof)
!print*,shape(sizescale)
!print*,shape(moist)

FuelFile = 'fuelfiles_RMSZ.dat'
grassfile = 'grass_DMSR.dat'

open(unit=1,file=FuelFile,form='unformatted',status='unknown')
write(1) Drhof,Dmoist,Dss,zheight
close(1)

open(unit=5,file=grassfile,form='unformatted',status='unknown')
write(5) gdepth,gmoisture,gss,grho
close(5)

if(windprofile.eq.0) then
  open(unit=3,file='windprofile.dat',form='unformatted',status='unknown')
  write(3) uavg,vavg,ustd,vstd
  close(3)
endif

if(iFIA.eq.1) then
  open(unit=4,file='FIA.dat',form='unformatted',status='unknown', access='stream')
  write(4) FIA
  close(4)
endif


open(unit=10,file='DUETInputs',form='formatted',status='unknown')
write(10,*) '&duetlist'
write(10,*) '                              '
write(10,*) '      inputprogram = ',inputprogram
write(10,*) '                              '
write(10,*) '!----------------------------!'
write(10,*) '! Grid variables'
write(10,*) '!----------------------------!'
write(10,*) '      nx = ',nx
write(10,*) '      ny = ',ny
write(10,*) '      nz = ',nz
write(10,*) '      dx = ',dx
write(10,*) '      dy = ',dy
write(10,*) '      zmax = ',zmax
write(10,*) '      PI = ',PI
write(10,*) '                               '
write(10,*) '!----------------------------!'
write(10,*) '! Fuel variables'
write(10,*) '!----------------------------!'
write(10,*) '      fueltotal = ',fueltotal
write(10,*) '      infuel = ',infuel
write(10,*) '      singlefuel = ',singlefuel
write(10,*) '      ntspecies = ',ntspecies
write(10,*) '      tfuelbins = ',tfuelbins
write(10,*) '      grassconstant = ',grassconstant
write(10,*) '      litterconstant = ',litterconstant
write(10,*) '      ngrass = ',duet_ngrass
write(10,*) '      grassstep = ',grassstep
write(10,*) '      gmoistoverride = ',gmoistoverride
write(10,*) '                                '
write(10,*) '      iFIA = ',iFIA
write(10,*) '      iFIAspecies = ',iFIAspecies
write(10,*) '                                '
write(10,*) '!----------------------------!'
write(10,*) '! Wind and Time variables'
write(10,*) '!----------------------------!'
write(10,*) '      windprofile = ',windprofile
write(10,*) '      randomwinds = ',randomwinds
write(10,*) '                              '
write(10,*) '      YearsSinceBurn = ',YearsSinceBurn
write(10,*) '      StepsPerYear = ',StepsPerYear
write(10,*) '                              '
write(10,*) '      relhum = ',relhum
write(10,*) '                              '
write(10,*) '!----------------------------!'
write(10,*) '! Import variables'
write(10,*) '!----------------------------!'
write(10,*) '      FuelFile = ',"'"//FuelFile//"'"
write(10,*) '      grassfile = ',"'"//grassfile//"'"
write(10,*) '                                '
write(10,*) '!----------------------------!'
write(10,*) '! Extra options'
write(10,*) '!----------------------------!'
write(10,*) '      litout = ',litout
write(10,*) '      controlseed = ',controlseed
write(10,*) '      seedchange = ',seedchange
write(10,*) '                              '
write(10,*) '/'
close(10)


end subroutine Duet_Inputs
