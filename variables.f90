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
!-----------------------------------------------------------------
module constant_variables
!-----------------------------------------------------------------
! Constant variables and arrays
!-----------------------------------------------------------------
implicit none
      
real  PI

end module constant_variables

module grid_variables
!-----------------------------------------------------------------
! Grid and topo variables and arrays
!-----------------------------------------------------------------
implicit none
      
integer nx,ny,nz
real    dx,dy,dz
real :: aa1=0.1
real :: ndatax=0.,ndatay=0.,datalocx=0.,datalocy=0. !JSM added ndatax, etc.
integer nfuel,zmax
real,allocatable:: rhof(:,:,:,:),sizescale(:,:,:,:),moist(:,:,:,:),fueldepth(:,:,:,:)
real,allocatable:: zs(:,:),zheight(:,:,:)
character:: topofile*50='flat'

end module grid_variables

module io_variables
!-----------------------------------------------------------------
! IO variables and arrays
!-----------------------------------------------------------------
implicit none
      
integer :: singlefuel=0,firetecshock=0
integer :: controlseed,n,seedchange
integer,allocatable :: seed(:)

end module io_variables

module infile_variables
!-----------------------------------------------------------------
! Variables for importing fuel files
!-----------------------------------------------------------------
implicit none
      
integer :: ifuelin=0
integer :: inx=0,iny=0,inz=0
real    :: idx=0.,idy=0.,idz=0.
real    :: iaa1=-1.
integer :: iintpr=0
integer :: infuel=0
real,allocatable:: irhof(:,:,:,:),iss(:,:,:,:),imoist(:,:,:,:),iafd(:,:,:,:)
real,allocatable:: izs(:,:),izheight(:,:,:)
character:: intopofile*50='flat' !JO
character:: rhoffile*50,moistfile*50,ssfile*50,afdfile*50

end module infile_variables

module baseline_variables
!-----------------------------------------------------------------
! Tree and groundfuel variables unique to baseline establishment
!-----------------------------------------------------------------
implicit none

integer:: igrass=0,itrees=0,ilitter=0
integer:: ngrass=1
real:: grassconstant=5.,litterconstant=5.,gmoistoverride=0.
real,allocatable:: grhof(:,:,:,:),gsizescale(:,:,:,:),gmoist(:,:,:,:),gfueldepth(:,:,:)
real,allocatable:: trhof(:,:,:,:),tsizescale(:,:,:,:),tmoist(:,:,:,:),tfueldepth(:,:,:)
real,allocatable:: lrhof(:,:,:,:),lsizescale(:,:,:,:),lmoist(:,:,:,:),lfueldepth(:,:,:)
character:: grassfile*50,treefile*50,litterfile*50

integer:: istem=0,ntspecies=1,tfuelbins=1,ntreefueltypes
real,allocatable:: tstemdensity(:),tlocation(:,:,:)
real,allocatable:: t1moisture(:,:),t1ss(:,:),t1bulkdensity(:,:)
real,allocatable:: t2moisture(:,:,:),t2ss(:,:,:),t2bulkdensity(:,:,:)
real,allocatable:: theight(:,:),tcrownbotheight(:,:),tcrownmaxheight(:,:),tcrowndiameter(:,:)
real,allocatable:: trhomicro(:),tdbh(:,:),tstemmoist(:),tbarkthick(:,:),tbarkmoist(:)
integer,allocatable:: ntrees(:),tspecies(:)
real,allocatable:: gdepth(:),grho(:),gss(:),gmoisture(:)
real,allocatable:: ldepth(:),lrho(:),lss(:),lmoisture(:)
integer,allocatable:: tdnx(:),tdny(:)

end module baseline_variables

module treatment_variables
!-----------------------------------------------------------------
! Types of treatments which can occur
!-----------------------------------------------------------------
implicit none
integer:: itreatment=0

!-----------------------------------------------------------------
! Slash Treatment variables
!-----------------------------------------------------------------

integer,allocatable:: sdnx(:),sdny(:)
real sdiameter,sheight
real sdepth
real sprho,smoist

end module treatment_variables

module duet_variables
!-----------------------------------------------------------------
! Surface fuel variables
!-----------------------------------------------------------------
implicit none

character :: speciesfile*100,winddatafile*100
integer :: windprofile=0,randomwinds=0
integer :: grassstep=1,periodTotal,litout
integer :: StepsPerYear=1,YearsSinceBurn=4
real :: relhum
real,allocatable:: vterminal(:),fuelSA(:),Froude(:),droptime(:)
real,allocatable:: leafdropfreq(:),decay(:),dragco(:),moistspec(:)
real,allocatable:: uavg(:),vavg(:),VAR(:,:),ustd(:),vstd(:)
real,allocatable:: ssspec(:),compact(:)
real,allocatable:: Umean(:,:,:),Vmean(:,:,:),Uvar(:,:,:),Vvar(:,:,:)
real,allocatable:: lrhofT(:,:,:,:),grhofT(:,:,:,:)
real,allocatable:: lafdT(:,:,:,:),gafdT(:,:,:,:)
real,allocatable:: lmoistT(:,:,:,:),gmoistT(:,:,:,:)
real,allocatable:: lssT(:,:,:,:),gssT(:,:,:,:)

end module duet_variables

module species_variables

!-----------------------------------------------------------------
! Variables for species database
!-----------------------------------------------------------------
implicit none

type :: read_species
  integer :: FIA_code
  character(len=8) :: sp_grp
  integer :: sp_grp_num
  character(len=30) :: species, genus, common_name
  real :: mass_avg,mass,surfarea_avg,surfarea
  integer :: dropperyear
  real :: decay,moist
  integer :: stepperyear
  real :: dragco,vterminal,froude,compact,sizescale
end type read_species

type :: read_species_grp
  real :: mass,surfarea
  integer :: dropperyear
  real :: decay,moist
  integer :: stepperyear
  real :: dragco,vterminal,froude,compact,sizescale
end type read_species_grp

type(read_species),dimension(290) :: SPECINFO
type(read_species_grp),dimension(10) :: SPECgroups
integer :: iFIA,iFIAspecies
integer,allocatable :: FIA(:)

end module species_variables



