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
integer nfuel,zmax
real,allocatable:: rhof(:,:,:,:),sizescale(:,:,:,:),moist(:,:,:,:),fueldepth(:,:,:,:)
real, allocatable :: trhof(:,:,:,:),tsizescale(:,:,:,:),tmoist(:,:,:,:),tfueldepth(:,:,:,:)
real,allocatable:: zheight(:,:,:)
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

module baseline_variables
!-----------------------------------------------------------------
! Tree and groundfuel variables unique to baseline establishment
!-----------------------------------------------------------------
implicit none

integer:: ngrass=1
integer :: ifuelin=0
integer :: infuel=0, fueltotal=0
real:: grassconstant=5.,litterconstant=5.,gmoistoverride=0.
real,allocatable:: grhof(:,:,:,:),gsizescale(:,:,:,:),gmoist(:,:,:,:),gfueldepth(:,:,:)
real,allocatable:: lrhof(:,:,:,:),lsizescale(:,:,:,:),lmoist(:,:,:,:),lfueldepth(:,:,:)
real,allocatable:: gdepth(:),grho(:),gss(:),gmoisture(:)
character:: grassfile*50,FuelFile*50

integer:: istem=0,ntspecies=1,tfuelbins=1,ntreefueltypes

end module baseline_variables


module duet_variables
!-----------------------------------------------------------------
! Surface fuel variables
!-----------------------------------------------------------------
implicit none

character :: speciesfile*100,winddatafile*100
integer :: windprofile=0,randomwinds=0
integer :: grassstep=1,periodTotal,litout
integer :: StepsPerYear=1,YearsSinceBurn=4,inputprogram=1
real :: relhum,densitythresh
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
real,allocatable :: FIAtemp(:)

end module species_variables

!-----------------------------------------------------------------
! Variables for FastFuels voxelated array inputs
!-----------------------------------------------------------------

module FF_variables

!integer, parameter :: int16 = selected_int_kind(16)

real :: winddirection,windvary
real,allocatable :: FFrhof(:,:,:),FFmoist(:,:,:)
real,allocatable :: surfrhof(:,:,:),surfdepth(:,:,:),surfmoist(:,:,:)
integer*2,allocatable :: FFspec(:,:,:)
integer*2,allocatable :: specarray(:)

end module FF_variables


