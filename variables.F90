!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! variables declares all the constant variables used throughout the 
! program
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
! Constant variables and arrays
!-----------------------------------------------------------------------
module constant_variables
  implicit none
      
  real :: PI = 3.14159 !YUM
  real :: tolerance = 1.e-8

end module constant_variables

!-----------------------------------------------------------------------
! Grid and topo variables and arrays
!-----------------------------------------------------------------------
module grid_variables
  implicit none
      
  integer :: nx,ny,nz
  real :: dx,dy,dz
  real :: aa1=0.1
  real :: ndatax=0.,ndatay=0.,datalocx=0.,datalocy=0. !JSM added ndatax, etc.
  integer :: nfuel
  real,allocatable :: rhof(:,:,:,:),sizescale(:,:,:,:)
  real,allocatable :: moist(:,:,:,:),fueldepth(:,:,:,:)
  real,allocatable :: zs(:,:),zheight(:,:,:)
  character :: topofile*50='flat'

end module grid_variables

!-----------------------------------------------------------------------
! IO variables and arrays
!-----------------------------------------------------------------------
module io_variables
  implicit none
        
  integer :: singlefuel=0,lreduced=0
  integer :: verbose=0,doubleprec=0
  integer :: controlseed,n,seedchange
  integer,allocatable :: seed(:)
  character:: workdir*255=''
  character :: filesep*1=''

end module io_variables

!-----------------------------------------------------------------------
! Variables for importing fuel files
!-----------------------------------------------------------------------
module infile_variables
  implicit none
        
  integer :: ifuelin=0
  integer :: inx=0,iny=0,inz=0
  integer :: iintpr=0
  integer :: infuel=0
  real :: idx=0.,idy=0.,idz=0.
  real :: iaa1=-1.
  real,allocatable :: irhof(:,:,:,:),iss(:,:,:,:)
  real,allocatable :: imoist(:,:,:,:),iafd(:,:,:,:)
  real,allocatable :: izs(:,:),izheight(:,:,:)
  character :: intopofile*50='flat' 
  character :: rhoffile*50,moistfile*50,ssfile*50,afdfile*50

end module infile_variables

!-----------------------------------------------------------------------
! Tree and groundfuel variables unique to baseline establishment
!-----------------------------------------------------------------------
module fuel_variables
  implicit none
  
  integer :: igrass=0,itrees=0,ilitter=0
  integer :: ngrass=0,duet_ngrass=0
  integer :: iFIA,iFIAspecies
  real :: grassconstant=5.,litterconstant=5.,gmoistoverride=0.
  real,allocatable :: grhof(:,:,:,:),gsizescale(:,:,:,:)
  real,allocatable :: gmoist(:,:,:,:),gfueldepth(:,:,:)
  real,allocatable :: trhof(:,:,:,:),tsizescale(:,:,:,:)
  real,allocatable :: tmoist(:,:,:,:),tfueldepth(:,:,:)
  real,allocatable :: lrhof(:,:,:,:),lsizescale(:,:,:,:)
  real,allocatable :: lmoist(:,:,:,:),lfueldepth(:,:,:)
  character:: treefile*250,newtreefile*250
  
  integer :: istem=0,ntspecies=0,tfuelbins=0,ntreefueltypes
  integer,allocatable :: ntrees(:),tspecies(:),tunique_species(:)
  integer,allocatable :: trhofmaxindex(:)
  real,allocatable :: tstemdensity(:),tlocation(:,:,:)
  real,allocatable :: t1moisture(:,:),t1ss(:,:),t1bulkdensity(:,:)
  real,allocatable :: t2moisture(:,:,:),t2ss(:,:,:),t2bulkdensity(:,:,:)
  real,allocatable :: theight(:,:),tcrownbotheight(:,:)
  real,allocatable :: tcrownmaxheight(:,:),tcrowndiameter(:,:)
  real,allocatable :: trhomicro(:),tdbh(:,:),tstemmoist(:)
  real,allocatable :: tbarkthick(:,:),tbarkmoist(:)
  real,allocatable :: trhofmax(:)
  real,allocatable :: gdepth(:),grho(:),gss(:),gmoisture(:)
  real,allocatable :: ldepth(:),lrho(:),lss(:),lmoisture(:)
  
  character :: command*50

end module fuel_variables
