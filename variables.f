!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! variables declares all the constant variables used throughout the 
! program
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      real    aa1,ndatax,ndatay,datalocx,datalocy !JSM added ndatax, etc.
      integer nfuel,zmax,singlefuel,ifiretecshock
      real,allocatable:: rhof(:,:,:,:),sizescale(:,:,:,:),moist(:,:,:,:),fueldepth(:,:,:,:)
      real,allocatable:: zs(:,:),zheight(:,:,:)
      character:: topofile*50
      
      end module grid_variables
      
      module infile_variables
      !-----------------------------------------------------------------
      ! Variables for importing fuel files
      !-----------------------------------------------------------------
      implicit none
            
      integer ifuelin
      integer inx,iny,inz
      real    idx,idy,idz
      real    iaa1
      integer iintpr
      integer infuel
      real,allocatable:: irhof(:,:,:,:),iss(:,:,:,:),imoist(:,:,:,:),iafd(:,:,:,:)
      real,allocatable:: izs(:,:),izheight(:,:,:)
      character:: rhoffile*50,moistfile*50,ssfile*50,afdfile*50
      
      end module infile_variables

      module baseline_variables
      !-----------------------------------------------------------------
      ! Tree and groundfuel variables unique to baseline establishment
      !-----------------------------------------------------------------
      implicit none
    
      integer:: igrass,itrees,ilitter 
      integer:: ngrass
      real:: grassconstant,litterconstant
      real,allocatable:: grhof(:,:,:,:),gsizescale(:,:,:,:),gmoist(:,:,:,:),gfueldepth(:,:,:)
      real,allocatable:: trhof(:,:,:,:),tsizescale(:,:,:,:),tmoist(:,:,:,:),tfueldepth(:,:,:)
      real,allocatable:: lrhof(:,:,:,:),lsizescale(:,:,:,:),lmoist(:,:,:,:),lfueldepth(:,:,:)
      character:: grassfile*50,treefile*50,litterfile*50

      integer:: istem,ntspecies,tfuelbins,ntreefueltypes
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
      integer:: itreatment
      
      !-----------------------------------------------------------------
      ! Slash Treatment variables
      !-----------------------------------------------------------------
      
      integer,allocatable:: sdnx(:),sdny(:)
      real sdiameter,sheight
      real sdepth
      real sprho,smoist

      end module treatment_variables
