!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! io contains the functions for reading input namelists and tree data 
! files and writing .dat files
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
subroutine fuellist_input
!-----------------------------------------------------------------
! namelist_input is a function reads in the user-defined variables 
! from the treelist then assigns those variables for use by the 
! trees softeware.
!-----------------------------------------------------------------
use grid_variables
use io_variables
use infile_variables
use fuels_create_variables
!use duet_variables, only : speciesfile,winddatafile,windprofile, &
!  grassstep,YearsSinceBurn,StepsPerYear,randomwinds,relhum, &
!  ustd,vstd,uavg,vavg,periodTotal,litout,densitythresh
!use species_variables
implicit none

!Local Variables

! Executable Code
open(unit=48,file= 'fuellist',form= 'formatted',status= 'old')

! Domain information
call QueryFuellist_integer('nx',nx,48,200)
call QueryFuellist_integer('ny',ny,48,200)
call QueryFuellist_integer('nz',nz,48,41)
call QueryFuellist_real('dx',dx,48,2.)
call QueryFuellist_real('dy',dy,48,2.)
call QueryFuellist_real('dz',dz,48,15.)
call QueryFuellist_real('aa1',aa1,48,0.1)
call QueryFuellist_integer('singlefuel',singlefuel,48,0)
call QueryFuellist_integer('lreduced',lreduced,48,0)
call QueryFuellist_string('topofile',topofile,48,'flat')
call QueryFuellist_integer('doubleprec',doubleprec,48,0)

! Data import from existing files
call QueryFuellist_integer('ifuelin',ifuelin,48,0)
if(ifuelin.eq.1) then
  call QueryFuellist_integer('inx',inx,48,nx)
  call QueryFuellist_integer('iny',iny,48,ny)
  call QueryFuellist_integer('inz',inz,48,nz)
  call QueryFuellist_real('idx',idx,48,dx)
  call QueryFuellist_real('idy',idy,48,dy)
  call QueryFuellist_real('idz',idz,48,dz)
  call QueryFuellist_real('iaa1',iaa1,48,aa1)
  call QueryFuellist_integer('infuel',infuel,48,1)
  call QueryFuellist_string('intopofile',intopofile,48,'flat')
  call QueryFuellist_string('rhoffile',rhoffile,48,'treesrhof.dat.orig')
  call QueryFuellist_string('ssfile',ssfile,48,'treesss.dat.orig')
  call QueryFuellist_string('moistfile',moistfile,48,'treesmoist.dat.orig')
  call QueryFuellist_string('afdfile',afdfile,48,'treesfueldepth.dat.orig')
endif
    
call QueryFuellist_integer('ilitter',ilitter,48,0)
call QueryFuellist_integer('igrass',igrass,48,0)
if(ilitter.eq.2) igrass = 1

! Grass

if(igrass.eq.1)then
  call QueryFuellist_integer('ngrass',ngrass,48,1)
  call QueryFuellist_real('grassconstant',grassconstant,48,5.)
  allocate(grho(ngrass))
  call QueryFuellist_real_array('grho',grho,ngrass,48,1.18)
  allocate(gmoisture(ngrass))
  call QueryFuellist_real_array('gmoisture',gmoisture,ngrass,48,0.06)
  allocate(gss(ngrass))
  call QueryFuellist_real_array('gss',gss,ngrass,48,0.0005)
  allocate(gdepth(ngrass))
  call QueryFuellist_real_array('gdepth',gdepth,ngrass,48,0.27)
endif

! Trees
call QueryFuellist_integer('itrees',itrees,48,0)
if(itrees.ge.1)then
  call QueryFuellist_integer('tfuelbins',tfuelbins,48,1)
  call QueryFuellist_string('treefile',treefile,48,'AJoseEglinTrees.txt')
  call QueryFuellist_real('ndatax',ndatax,48,nx*dx)
  call QueryFuellist_real('ndatay',ndatay,48,ny*dy)
  call QueryFuellist_real('datalocx',datalocx,48,0.)
  call QueryFuellist_real('datalocy',datalocy,48,0.)
  call find_numspecies
  allocate(trhofmax(ntspecies))
  call QueryFuellist_real_array('trhofmax',trhofmax,ntspecies,48,0.)
endif

! Litter

if(ilitter.eq.1) then
  call QueryFuellist_real('litterconstant',litterconstant,48,0.0)
  allocate(lrho(ntspecies*tfuelbins))
  call QueryFuellist_real_array('lrho',lrho,ntspecies*tfuelbins,48,1.18)
  allocate(lmoisture(ntspecies*tfuelbins))
  call QueryFuellist_real_array('lmoisture',lmoisture,ntspecies*tfuelbins,48,0.06)
  allocate(lss(ntspecies*tfuelbins))
  call QueryFuellist_real_array('lss',lss,ntspecies*tfuelbins,48,0.0005)
  allocate(ldepth(ntspecies*tfuelbins))
  call QueryFuellist_real_array('ldepth',ldepth,ntspecies*tfuelbins,48,0.05)
!elseif(ilitter.eq.2)then  ! DUET
!  call QueryFuellist_integer('windprofile',windprofile,48,2)
!  call QueryFuellist_integer('YearsSinceBurn',YearsSinceBurn,48,1)
!  call QueryFuellist_integer('StepsPerYear',StepsPerYear,48,1)
!  call QueryFuellist_real('relhum',relhum,48,0.1)
!  call QueryFuellist_integer('grassstep',grassstep,48,1)
!  call QueryFuellist_integer('iFIA',iFIA,48,1)
!  call QueryFuellist_integer('iFIAspecies',iFIAspecies,48,1)
!  call QueryFuellist_integer('litout',litout,48,0)
!  call QueryFuellist_real('gmoistoverride',gmoistoverride,48,0.0)
!  call QueryFuellist_real('densitythresh',densitythresh,48,0.5)
!  if(windprofile.eq.0)then
!    periodTotal=YearsSinceBurn*StepsPerYear
!    allocate(uavg(periodTotal))
!    call QueryFuellist_real_array('uavg',uavg,periodTotal,48,0.)
!    allocate(vavg(periodTotal))
!    call QueryFuellist_real_array('vavg',vavg,periodTotal,48,0.)
!    allocate(ustd(periodTotal))
!    call QueryFuellist_real_array('ustd',ustd,periodTotal,48,0.)
!    allocate(vstd(periodTotal))
!    call QueryFuellist_real_array('vstd',vstd,periodTotal,48,0.)
!  elseif(windprofile.eq.2)then
!    call QueryFuellist_integer('randomwinds',randomwinds,48,2)
!  endif
!  if(iFIA.eq.0)then
!    call QueryFuellist_string('speciesfile',speciesfile,48,'speciesfile.dat')
!  elseif(iFIA.eq.1)then 
!    allocate(FIA(infuel+ntspecies*tfuelbins))
!    call QueryFuellist_integer_array('FIA',FIA,infuel+ntspecies*tfuelbins,48,1)
!    call define_species_variables
!  endif
endif

! Verbose
call QueryFuellist_integer('verbose',verbose,48,0)
     
close(48)

end subroutine fuellist_input

subroutine output_fuel
!-----------------------------------------------------------------
! output is a function which writes the .dat files for use in 
! FIRETEC or QUIC-Fire
!-----------------------------------------------------------------
use grid_variables
use io_variables

implicit none

! Local Variables
integer ift,i,j,k,lfuel,lz
integer,dimension(nfuel):: nonzero
real(8),allocatable :: rho(:,:,:,:),h20(:,:,:,:)
real(8),allocatable :: sss(:,:,:,:),afd(:,:,:,:)

! Executable Code
if (singlefuel.eq.1) then
  do i=1,nx
    do j=1,ny
      do k=1,nz
        if(sum(rhof(:,i,j,k)).gt.0)then
          sizescale(1,i,j,k)=sum(sizescale(:,i,j,k)*rhof(:,i,j,k))/sum(rhof(:,i,j,k))
          moist(1,i,j,k)=sum(moist(:,i,j,k)*rhof(:,i,j,k))/sum(rhof(:,i,j,k))
          fueldepth(1,i,j,k)=sum(fueldepth(:,i,j,k)*rhof(:,i,j,k))/sum(rhof(:,i,j,k))
          rhof(1,i,j,k)=sum(rhof(:,i,j,k))
          if(nfuel.gt.1) rhof(2:nfuel,i,j,k)=0.
        endif
      enddo
    enddo
  enddo
endif
  
nonzero(:) = 1
lfuel = 1
do ift=1,nfuel
  if (sum(rhof(ift,:,:,:)).le.0.0001*sum(rhof)) then
    nonzero(ift) = 0
  else
    do k=1,nz
      if (sum(rhof(ift,:,:,k)).gt.0) lfuel = max(lfuel,k)
    enddo
  endif
enddo
if (lreduced.eq.1) then
  lz=lfuel
else
  lz=nz
endif

if(doubleprec.eq.1) then
  print*,'Converting to Double Precision'
  allocate(rho(nfuel,nx,ny,nz))
  allocate(h20(nfuel,nx,ny,nz))
  allocate(sss(nfuel,nx,ny,nz))
  allocate(afd(nfuel,nx,ny,nz))

  rho = 0.
  h20 = 0.
  sss = 0.
  afd = 0.

  rho = dble(rhof)
  h20 = dble(moist)
  sss = dble(sizescale)
  afd = dble(fueldepth)

  print*,'Exporting data to .dat files'
  call write_file_dbl('treesrhof.dat',rho(:,:,:,1:lz),nfuel,nx,ny,lz,nonzero)
  call write_file_dbl('treesmoist.dat',h20(:,:,:,1:lz),nfuel,nx,ny,lz,nonzero)
  call write_file_dbl('treesss.dat',sss(:,:,:,1:lz),nfuel,nx,ny,lz,nonzero)
  call write_file_dbl('treesfueldepth.dat',afd(:,:,:,1:lz),nfuel,nx,ny,lz,nonzero)
  
else

  print*,'Exporting data to .dat files'
  call write_file('treesrhof.dat',rhof(:,:,:,1:lz),nfuel,nx,ny,lz,nonzero)
  call write_file('treesmoist.dat',moist(:,:,:,1:lz),nfuel,nx,ny,lz,nonzero)
  call write_file('treesss.dat',sizescale(:,:,:,1:lz),nfuel,nx,ny,lz,nonzero)
  call write_file('treesfueldepth.dat',fueldepth(:,:,:,1:lz),nfuel,nx,ny,lz,nonzero)

endif

print*,'Your nfuel is',int(sum(nonzero(:)))
print*,'Your lfuel is',lfuel

end subroutine output_fuel

subroutine write_file(filename,dataname,nfuel,nx,ny,nz,nonzero)
!-----------------------------------------------------------------
! write_file writes data to specified file
!-----------------------------------------------------------------
implicit none

! Local Variables
character(len=*),intent(in) :: filename
integer,intent(in) :: nfuel,nx,ny,nz
integer,intent(in) :: nonzero(nfuel)
real,intent(in) :: dataname(nfuel,nx,ny,nz)

integer :: ift

! Executable Code
open (1,file= filename,form= 'unformatted',status= 'unknown')
do ift=1,nfuel
  if (nonzero(ift).ne.0)  write (1) dataname(ift,:,:,:)
enddo
close (1)

end subroutine write_file

subroutine write_file_dbl(filename,dataname,nfuel,nx,ny,nz,nonzero)
  !-----------------------------------------------------------------
  ! write_file writes data to specified file
  !-----------------------------------------------------------------
  implicit none
  
  ! Local Variables
  character(len=*),intent(in) :: filename
  integer,intent(in) :: nfuel,nx,ny,nz
  integer,intent(in) :: nonzero(nfuel)
  real(8),intent(in) :: dataname(nfuel,nx,ny,nz)
  
  integer :: ift
  
  ! Executable Code
  open (1,file= filename,form= 'unformatted',status= 'unknown')
  do ift=1,nfuel
    if (nonzero(ift).ne.0)  write (1) dataname(ift,:,:,:)
  enddo
  close (1)
  
  end subroutine write_file_dbl
  
!-------------------------------------------
! Write New Fuellist based on old params
!-------------------------------------------
subroutine output_fuellist

use grid_variables
use io_variables
use infile_variables
use fuels_create_variables
!use duet_variables
!use species_variables

!stuff for writing new treelist
character(len=100) :: newtreelistname
integer :: ii

open (2222,file= 'fuellist_'//newtreefile,form= 'formatted',status= 'unknown')
newtreelistname = newtreefile//'_treelist.txt'
write(2222,'(A10)') '&fuellist'
write(2222,'(A35)') '! ----------------------------------'
write(2222,'(A35)') '! FIRETEC domain info'
write(2222,'(A35)') '! ----------------------------------'
write(2222,'(A4,I6.2)')     'nx = ',nx
write(2222,'(A4,I6.2)')     'ny = ',ny
write(2222,'(A4,I6.2,A50)') 'nz = ',nz,'! Size of HIGRAD/FIRETEC grid [cells]'
write(2222,'(A4,F6.2)')     'dx = ',dx
write(2222,'(A4,F6.2)')     'dy = ',dy
write(2222,'(A4,F6.2,A50)') 'dz = ',dz,'! Grid Resolution [m]'
write(2222,'(A6,F6.2,A50)') 'aa1 = ',aa1,'! Vertical stretching component [default=0.1]'
write(2222,'(A15,I2.1,A60)') 'singlefuel = ',singlefuel,'! Flag forcing single fuel type instead of multiple fuels'
write(2222,'(3A30)') "topofile = '"//TRIM(topofile)//"' 	     "

if(ifuelin.eq.1) then
    write(2222,'(A35)')'! ----------------------------------'
    write(2222,'(A35)')'! Data import from existing files'
    write(2222,'(A35)')'! ----------------------------------'
    write(2222,'(A25)')'ifuelin = 0 ! Fuel read-in flag '
    write(2222,'(A6)')'inx = 0'
    write(2222,'(A6)')'iny = 0'
    write(2222,'(A30)')'inz = 0 ! Size of HIGRAD/FIRETEC grid'
    write(2222,'(A6)')'idx = 0'
    write(2222,'(A6)')'idy = 0'
    write(2222,'(A30)')'idz = 0 ! Grid Resolution [m]'
    write(2222,'(A50)')'iaa1 = 0.1 ! Vertical stretching component'
    write(2222,'(A30)')'infuel = 1 ! Number of Fuel Types'
    write(2222,'(A15)')'intopofile =  " " '
    write(2222,'(A15)')'rhoffile =  " " '
    write(2222,'(A15)')'ssfile   =  " " '
    write(2222,'(A15)')'moistfile =  " " '
    write(2222,'(A15)')'afdfile  =  " " '
endif

write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A35)')'! Input trees dataset info'
write(2222,'(A35)')'! ----------------------------------'
if (itrees.eq.0) then 
    write(2222,'(A25)')'itrees = 0 ! Trees flag '
else
    write(2222,'(A25)')'itrees = 2 ! Trees flag '
endif

if(itrees.ge.1)then
    write(2222,'(A12,I1)') "tfuelbins = ",tfuelbins
    write(2222,*) "treefile = '",TRIM(newtreelistname)//'_treelist.txt',"'"
    !write(2222,'(A10,A30,A1)') "treefile= '",TRIM(newtreelistname),"'"
    write(2222,'(A10,F6.1)') "ndatax = ",nx*dx
    write(2222,'(A10,F6.1)') "ndatay = ",ny*dy
    write(2222,'(A10,F6.1)') "datalocx = ",datalocx
    write(2222,'(A10,F6.1)') "datalocy = ",datalocy
endif 

write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A35)')'! Litter switch'
write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A15,I1)')'ilitter = ',ilitter
if(ilitter.eq.1) then
    write(2222,'(A35)')'!ilitter eq 1 (BASIC) info'
    write(2222,'(A25,F4.2)')'litterconstant = ',litterconstant
    write(2222,'(A10)', Advance = 'No' ) 'lrho = '
    do ii = 1, ntspecies*tfuelbins
      write(2222, '( F5.3, A1 )', Advance = 'No' ) lrho(ii)
      if (ii.ne.ntspecies*tfuelbins) write(2222,'(A1)', Advance = 'No') ','
    end do
    write(2222,'(/,A15)', Advance = 'No' ) 'lmoisture = ' 
    do ii = 1, ntspecies*tfuelbins
      write(2222, '( F4.3, A1 )', Advance = 'No' ) lmoisture(ii)
      if (ii.ne.ntspecies*tfuelbins) write(2222,'(A1)', Advance = 'No') ','
    end do
    write(2222,'(/,A6)', Advance = 'No' ) 'lss = ' 
    do ii = 1, ntspecies*tfuelbins
      write(2222, '( F6.5, A1 )', Advance = 'No' ) lss(ii)
      if (ii.ne.ntspecies*tfuelbins) write(2222,'(A1)', Advance = 'No') ','
    end do
    write(2222,'(/,A10)', Advance = 'No' ) 'ldepth = ' 
    do ii = 1, ntspecies*tfuelbins
      write(2222, '( F4.3, A1 )', Advance = 'No' ) ldepth(ii)
      if (ii.ne.ntspecies*tfuelbins) write(2222,'(A1)', Advance = 'No') ','
    end do
    
!elseif(ilitter.eq.2)then  ! DUET
!    write(2222,'(A35)')'!ilitter eq 2 (DUET) info  '
!    write(2222,'(A25,I1)')'windprofile = ',windprofile
!    write(2222,'(A25,I4.1)')'YearsSinceBurn = ',YearsSinceBurn
!    write(2222,'(A25,I4.1)')'StepsPerYear = ',StepsPerYear
!    write(2222,'(A25,F6.2)')'relhum = ',relhum
!    write(2222,'(A25,I3.1)')'grassstep = ',grassstep
!    write(2222,'(A25,I3.1)')'iFIA = ',iFIA
!    write(2222,'(A25,I3.1)')'iFIAspecies = ',iFIAspecies
!    write(2222,'(A25,I4.1)')'litout = ',litout
!    write(2222,'(A25,F3.1)') 'gmoistoverride = ',gmoistoverride
!    if(windprofile.eq.0)then
!        periodTotal=YearsSinceBurn*StepsPerYear
!        write(2222,'(A10)', Advance = 'No' ) 'uavg = '
!        do ii = 1, periodTotal
!          write(2222, '( F4.3, A1 )', Advance = 'No' ) uavg(ii) 
!          if (ii.ne.periodTotal) write(2222,'(A1)', Advance = 'No') ','
!        end do
!            write(2222,'(/,A10)', Advance = 'No' ) 'vavg = '
!        do ii = 1, periodTotal
!          write(2222, '( F4.3, A1 )', Advance = 'No' ) vavg(ii) 
!          if (ii.ne.periodTotal) write(2222,'(A1)', Advance = 'No') ','
!        end do
!            write(2222,'(/,A10)', Advance = 'No' ) 'ustd = '
!        do ii = 1, periodTotal
!          write(2222, '( F4.3, A1 )', Advance = 'No' ) ustd(ii)  
!          if (ii.ne.periodTotal) write(2222,'(A1)', Advance = 'No') ','
!        end do
!            write(2222,'(/,A10)', Advance = 'No' ) 'vstd = '
!        do ii = 1, periodTotal
!          write(2222, '( F4.3, A1 )', Advance = 'No' ) vstd(ii)  
!          if (ii.ne.periodTotal) write(2222,'(A1)', Advance = 'No') ','
!        end do
!    elseif(windprofile.eq.2)then
!        write(2222,'(A25,I2.1)')'randomwinds = ',randomwinds
!    endif
!
!    if (iFIA.eq.0)then
!        write(2222,*) "speciesfile = '",TRIM(speciesfile),"'"
!    elseif(iFIA.eq.1)then
!        write(2222,'(A5)', Advance = 'No' )'FIA = '
!        do ii = 1, ntspecies
!          write(2222, '( I4.3, A1 )', Advance = 'No' ) FIA(ii)
!          if (ii.ne.ntspecies) write(2222,'(A1)', Advance = 'No') ','
!        end do
!    endif
endif

write(2222,'(/,A35)')'! ----------------------------------'
write(2222,'(A35)')'! Grass switch'
write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A25,I1)')'igrass = ',igrass
if(igrass.eq.1)then
    write(2222,'(A35)')'!igrass options'
    write(2222,'(A25,I2.1)')'ngrass = ',ngrass
    write(2222,'(A25,F4.2)')'grassconstant = ',grassconstant
    write(2222,'(A10)', Advance = 'No' ) 'grho = '
    do ii = 1, ngrass
      write(2222, '( F8.4, A1 )', Advance = 'No' ) grho(ii) 
      if (ii.ne.ngrass) write(2222,'(A1)', Advance = 'No') ','
    end do
    write(2222,'(/,A12)', Advance = 'No' ) 'gmoisture = '
    do ii = 1, ngrass
      write(2222, '( F4.3, A1 )', Advance = 'No' ) gmoisture(ii) 
      if (ii.ne.ngrass) write(2222,'(A1)', Advance = 'No') ','
    end do
    write(2222,'(/,A10)', Advance = 'No' ) 'gss = '
    do ii = 1, ngrass
      write(2222, '( F6.4, A1 )', Advance = 'No' ) gss(ii)  
      if (ii.ne.ngrass) write(2222,'(A1)', Advance = 'No') ','
    end do
    write(2222,'(/,A10)', Advance = 'No' ) 'gdepth = '
    do ii = 1, ngrass
      write(2222, '( F4.3, A1 )', Advance = 'No' ) gdepth(ii) 
      if (ii.ne.ngrass) write(2222,'(A1)', Advance = 'No') ','
    end do
endif

write(2222,'(/,A35)')'! ----------------------------------'
write(2222,'(A50)')'! Option to output tree list and fuellist'
write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A75)')'verbose = 0 !flag to output treelist and fuellist from resulting run: 0=No ; 1=Yes'

close(2222)
  

end subroutine output_fuellist
