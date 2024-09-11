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

! This is a list of subroutines for the DUET program.
! Author:  Jenna Sjunneson McDanold 7/2024

module DUETio

  type :: domaininfo
    integer :: nx,ny,nz,ns,nt,ysb,spy,ng
    real :: dx,dy,dz
  end type domaininfo

  type(domaininfo) :: domain

  type :: inarrays
    real,allocatable :: trhof(:,:,:,:),zheight(:,:,:)
    real,allocatable :: moist(:,:,:,:),depth(:,:,:)
    real,allocatable :: sscale(:,:,:,:)
    character*50 :: speciesfile  
  end type inarrays
    
  type(inarrays) :: inarray

  type :: usearrays
    real,allocatable :: lrho(:,:,:,:),lsss(:,:,:,:),lafd(:,:,:,:),lh20(:,:,:,:)
  end type usearrays

  type(usearrays) :: litter
  
  type :: outarrays
    real,allocatable :: srho(:,:,:),ssss(:,:,:),safd(:,:,:),sh20(:,:,:) ! surface arrays (litter and grass)
    real,allocatable :: lrho(:,:,:),lsss(:,:,:),lafd(:,:,:),lh20(:,:,:) ! litter arrays
    real,allocatable :: frho(:,:,:,:),fsss(:,:,:,:),fafd(:,:,:),fh20(:,:,:,:) ! full arrays (litter and grass and trees)
  end type outarrays
  
  type(outarrays) :: outarray
  
  type :: specdetails
    real,allocatable :: comp(:),decy(:),frde(:),fuSA(:),drag(:)
    real,allocatable :: grho(:),vter(:),ssss(:),fh20(:),dept(:),fuMA(:)
    integer,allocatable :: drop(:),step(:)
  end type specdetails

  type(specdetails) :: species, grasses

  type :: windinfo
    real,allocatable :: uavg(:),vavg(:),uvar(:),vvar(:)
  end type windinfo

  type(windinfo) :: windarray

  type :: duetvariables
    integer :: windprofile=2,winddirection=0,windvary=90,randomwinds=3
    integer :: iFIA=1,iFIAspecies=1,litout=1
    integer :: controlseed=1
    real :: relhum=0.1,gmoistoverride=0,densitythresh=1.5
    real :: grassconstant=5.0,litterconstant=5.0,seedchange
  end type duetvariables

  type(duetvariables) :: duetvars

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
  
  integer,allocatable :: specarray(:)
  
end module DUETio
    
!---------------------------------------------------------------------!

module fillio

contains 
  subroutine TR_filldomain

    use DUETio

    character*50 :: speciesfile,filepath
    logical :: exists

    print*,'DUET begun...'

    inquire(file='../Inputs/fuellist', exist=exists)
    if(exists) then
      filepath = '../Inputs/fuellist'
    else
      inquire(file='fuellist', exist=exists)
      if(exists) then
        filepath = 'fuellist'
      endif
    endif

    open(unit=48,file=filepath,form= 'formatted',status= 'old')

    call QueryFuellist_real('dx',domain%dx,48,2.0)
    call QueryFuellist_real('dy',domain%dy,48,2.0)
    call QueryFuellist_real('dz',domain%dz,48,1.0)

    call QueryFuellist_integer('YearsSinceBurn',domain%ysb,48,1)
    call QueryFuellist_integer('StepsPerYear',domain%spy,48,1)
    call QueryFuellist_integer('ngrass',domain%ng,48,1)

    domain%nt = domain%ysb * domain%spy

    call alloc_init

    call QueryFuellist_integer('windprofile',duetvars%windprofile,48,2)

    call QueryFuellist_real('relhum',duetvars%relhum,48,0.1)

    call QueryFuellist_real('grassconstant',duetvars%grassconstant,48,5.0)
    call QueryFuellist_real('litterconstant',duetvars%litterconstant,48,5.0)
    call QueryFuellist_integer_array('grassstep',grasses%step,domain%ng,48,1)
    call QueryFuellist_real_array('grho',grasses%grho,domain%ng,48,1.1766)
    call QueryFuellist_real_array('gmoisture',grasses%fh20,domain%ng,48,0.06)
    call QueryFuellist_real_array('gss',grasses%ssss,domain%ng,48,0.0005)
    call QueryFuellist_real_array('gdepth',grasses%dept,domain%ng,48,0.27)

    call QueryFuellist_integer('iFIA',duetvars%iFIA,48,1)
    call QueryFuellist_integer('iFIAspecies',duetvars%iFIAspecies,48,1)

    call QueryFuellist_integer('litout',duetvars%litout,48,0)
    call QueryFuellist_real('gmoistoverride',duetvars%gmoistoverride,48,0.0)
    call QueryFuellist_real('densitythresh',duetvars%densitythresh,48,0.5)
    
    if(windprofile.eq.0)then
      call QueryFuellist_real_array('uavg',windarray%uavg,domain%nt,48,0.)
      call QueryFuellist_real_array('vavg',windarray%vavg,domain%nt,48,0.)
      call QueryFuellist_real_array('ustd',windarray%uvar,domain%nt,48,0.)
      call QueryFuellist_real_array('vstd',windarray%vvar,domain%nt,48,0.)
    elseif(windprofile.eq.2)then
      call QueryFuellist_integer('randomwinds',duetvars%randomwinds,48,3)
    endif

    if(iFIA.eq.0)then
      call QueryFuellist_string('speciesfile',speciesfile,48,'speciesfile.dat')
      call usespeciesfile(speciesfile)
    elseif(iFIA.eq.1)then 
      call QueryFuellist_integer_array('FIA',specarray,domain%ns,48,1)
    endif

    close(48)

  end subroutine TR_filldomain

  !---------------------------------------------------------------------!

  subroutine alloc_init

    use DUETio

    integer :: nx,ny,nt,ng,YSB

    nz = domain%nz
    nx = domain%nx
    ny = domain%ny
    nt = domain%nt
    ng = domain%ng
    ns = domain%ns
    YSB = domain%ysb

    allocate(litter%lrho(ns,nx,ny,nt))
    allocate(litter%lsss(ns,nx,ny,nt))
    allocate(litter%lafd(ns,nx,ny,nt))
    allocate(litter%lh20(ns,nx,ny,nt))

    allocate(outarray%srho(ns+ng,nx,ny))
    allocate(outarray%ssss(ns+ng,nx,ny))
    allocate(outarray%safd(ns+ng,nx,ny))
    allocate(outarray%sh20(ns+ng,nx,ny))

    allocate(outarray%lrho(ns,nx,ny))
    allocate(outarray%lsss(ns,nx,ny))
    allocate(outarray%lafd(ns,nx,ny))
    allocate(outarray%lh20(ns,nx,ny))

    allocate(outarray%frho(2*ns+ng,nx,ny,nz))
    allocate(outarray%fsss(2*ns+ng,nx,ny,nz))
    allocate(outarray%fafd(2*ns+ng,nx,ny))
    allocate(outarray%fh20(2*ns+ng,nx,ny,nz))

    outarray%srho = 0.0 ! surface fuels = grass and litter
    outarray%ssss = 0.0 ! surface fuels = grass and litter
    outarray%safd = 0.0 ! surface fuels = grass and litter
    outarray%sh20 = 0.0 ! surface fuels = grass and litter

    outarray%lrho = 0.0 ! just litter
    outarray%lsss = 0.0 ! just litter
    outarray%lafd = 0.0 ! just litter
    outarray%lh20 = 0.0 ! just litter

    outarray%frho = 0.0 ! full arrays = grass, litter, and trees
    outarray%fsss = 0.0 ! full arrays = grass, litter, and trees
    outarray%fafd = 0.0 ! full arrays = grass, litter, and trees
    outarray%fh20 = 0.0 ! full arrays = grass, litter, and trees

    allocate(species%comp(ns))
    allocate(species%decy(ns))
    allocate(species%frde(ns))
    allocate(species%fuSA(ns))
    allocate(species%vter(ns))
    allocate(species%ssss(ns))
    allocate(species%fh20(ns))
    allocate(species%drop(ns))
    allocate(species%step(ns))
    allocate(species%drag(ns))
    allocate(species%fuMA(ns))

    species%comp = 0.0
    species%decy = 0.0
    species%frde = 0.0
    species%fuSA = 0.0
    species%vter = 0.0
    species%ssss = 0.0
    species%fh20 = 0.0
    species%drop = 0.0
    species%step = 0.0
    species%drag = 0.0
    species%fuMA = 0.0

    allocate(grasses%decy(ng))
    allocate(grasses%fh20(ng))
    allocate(grasses%step(ng))
    allocate(grasses%grho(ng))
    allocate(grasses%dept(ng))
    allocate(grasses%ssss(ng))

    grasses%decy  = 1.0
    grasses%fh20  = 0.1
    grasses%step  = 1
    grasses%grho  = 0.0
    grasses%dept  = 0.1
    grasses%ssss  = 0.0005
    
    allocate(windarray%uavg(YSB))
    allocate(windarray%vavg(YSB))
    allocate(windarray%uvar(YSB))
    allocate(windarray%vvar(YSB))
    
    windarray%uavg = 0.0
    windarray%vavg = 0.0
    windarray%uvar = 0.0
    windarray%vvar = 0.0

    allocate(specarray(ns))

    call makeFIAfile
    
  end subroutine alloc_init

  !---------------------------------------------------------------------!

  subroutine makeFIAfile
    use DUETio, only : SPECINFO,SPECgroups

    integer :: i,j,ct,fdrop=0,fstep=0
    real :: fsa=0,fdecay=0,fmoist=0,fss=0
    real :: fdrag=0,fvterm=0,ffrou=0,fcomp=0

    open(unit=5, file='FIA_FastFuels_fin_fulllist_populated.txt', &
      status='old', access='sequential', form='formatted')
      do i=1,290
        read(5,*) SPECINFO(i)
      enddo
    close(5)

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
    
  end subroutine makeFIAfile

  !---------------------------------------------------------------------!

  subroutine usespeciesfile(speciesfile)

    use DUETio, only : domain,species

    character*50, intent(in) :: speciesfile

    real :: dragslope,dragb
    integer :: s

    print*,'!----!----!----!----!----!----!----!----!----!----!----!----!----!----!'
    print*,'YOU ARE NOT USING THE FIA DATABASE - SPECIES FILE PROVIDED SEPARATELY'
    print*,'!----!----!----!----!----!----!----!----!----!----!----!----!----!----!'
    
    dragslope = 1.4/(sqrt(70.0) - 1.0)
    dragb = 0.6 - dragslope
    
    open (99,file=speciesfile)
    do s=1,domain%ns
      read(99,*) species%fuMA(s),species%fuSA(s),species%drop(s),species%decy(s),species%step(s)
      species%fuMA(s)=species%fuMA(s)/1000.
      species%fuSA(s)=species%fuSA(s)/10000.
      species%drag(s)=dragslope*species%fuSA(s) + dragb !JSM
      ! Terminal velocity
      species%vter(s)=sqrt((2.*species%fuMA(s)*9.81)/(species%drag(s)*1.225*species%fuSA(s)))
      ! 9.81 is acceleration of gravity
      ! 1.225 is density of air at STP
      species%frde(s) = 0.01/(sqrt(9.81*sqrt(species%fuSA(s))))
      species%fh20(s)=100
      species%comp(s) = 0.5
    enddo
    close(99)

  end subroutine usespeciesfile

end module fillio




    