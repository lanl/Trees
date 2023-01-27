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
subroutine namelist_input
!-----------------------------------------------------------------
! namelist_input is a function reads in the user-defined variables 
! from the treelist then assigns those variables for use by the 
! trees softeware.
!-----------------------------------------------------------------
use grid_variables
use io_variables
use infile_variables
use baseline_variables
use treatment_variables
use duet_variables, only : speciesfile,winddatafile,windprofile, &
  grassstep,YearsSinceBurn,StepsPerYear,randomwinds,relhum, &
  ustd,vstd,uavg,vavg,periodTotal,litout
use species_variables
implicit none

!Local Variables
namelist/fuellist/ &
   nx,ny,nz,dx,dy,dz,aa1,topofile, &
   singlefuel,firetecshock,intopofile, &
   ifuelin,rhoffile,moistfile,ssfile,afdfile, &
   inx,iny,inz,idx,idy,idz,iaa1,infuel, &
   igrass,ngrass,grassconstant,grassfile, &
   itrees,ntspecies,iFIA,iFIAspecies, &
   tfuelbins,tdnx,tdny,treefile,istem, &
   ndatax,ndatay,datalocx,datalocy, & !JSM added for populate function
   ilitter,litterconstant,litterfile, &
   speciesfile,randomwinds,relhum,litout, &
   controlseed,seedchange, &
   winddatafile,windprofile,grassstep, &
   YearsSinceBurn,StepsPerYear, &
   itreatment,sdnx,sdny, &
   sdiameter,sheight,sprho,smoist,sdepth

 namelist/speciesdata/ &
   FIA

 namelist/litterdata/ &
   lrho,lmoisture,lss,ldepth

 namelist/grassdata/ &
   grho,gmoisture,gss,gdepth,gmoistoverride

 namelist/winddata/ &
   uavg,vavg,ustd,vstd


! Executable Code
! Area of interest arrays need to be allocated before calling namelist
allocate(tdnx(2)); tdnx(:)=0 ! Array of the cell range (x)  where the trees are applied
allocate(tdny(2)); tdny(:)=0 ! Array of the cell range (x)  where the trees are applied
allocate(sdnx(2)); sdnx(:)=0 ! Array of the cell range (x)  where the treatment is applied
allocate(sdny(2)); sdny(:)=0 ! Array of the cell range (x)  where the treatment is applied

!if(iFIA.eq.1.and.itrees.ne.7) call find_numspecies
!allocate(FIA(ntspecies))

open(unit=15,file='fuellist',form='formatted',status='old')
     read (15,nml=fuellist)
     if (ntspecies.eq.0) then
        iFIA = 0
        allocate(FIA(1*tfuelbins))        
     else
      allocate(FIA(ntspecies*tfuelbins))
      read (15,nml=speciesdata)
    endif
    
     if (ilitter.gt.0.and.ilitter.ne.2) then
       allocate(lrho(ntspecies*tfuelbins))
       allocate(lmoisture(ntspecies*tfuelbins))
       allocate(lss(ntspecies*tfuelbins))
       allocate(ldepth(ntspecies*tfuelbins))

       read (15,nml=litterdata)
     endif

     if (igrass.eq.1) then
       allocate(grho(ngrass))
       allocate(gmoisture(ngrass))
       allocate(gss(ngrass))
       allocate(gdepth(ngrass))

       read (15,nml=grassdata)

     endif

     if (windprofile.eq.0) then
       periodTotal=YearsSinceBurn*StepsPerYear
       allocate(uavg(periodTotal))
       allocate(vavg(periodTotal))
       allocate(ustd(periodTotal))
       allocate(vstd(periodTotal))

       read (15,nml=winddata)
     
     endif

close(15)


if (controlseed.eq.1) then
  call random_seed(size=n)
  allocate(seed(n))
  seed=seedchange
  call random_seed(put=seed)
  deallocate(seed)
endif

! Corrections for if variables not specified on namelist
if (tdnx(1).eq.0) then
  tdnx(1) = 1
  tdnx(2) = dx*nx
  tdny(1) = 1
  tdny(2) = dy*ny
endif
if (sdnx(1).eq.0) then
  sdnx(1) = 1
  sdnx(2) = dx*nx
  sdny(1) = 1
  sdny(2) = dy*ny
endif
tdnx(1)=floor(tdnx(1)/dx+1)
tdnx(2)=ceiling(tdnx(2)/dx)
tdny(1)=floor(tdny(1)/dy+1)
tdny(2)=ceiling(tdny(2)/dy)
sdnx(1)=floor(sdnx(1)/dx+1)
sdnx(2)=ceiling(sdnx(2)/dx)
sdny(1)=floor(sdny(1)/dy+1)
sdny(2)=ceiling(sdny(2)/dy)

!JO

if (iFIA.eq.1) then
  call define_species_variables
endif


if (ndatax.eq.0) ndatax=nx*dx
if (ndatay.eq.0) ndatay=ny*dy

if (ifuelin.eq.1) then
  if(inx.eq.0) inx=nx
  if(iny.eq.0) iny=ny
  if(inz.eq.0) inz=nz
  if(idx.eq.0) idx=dx
  if(idy.eq.0) idy=dy
  if(idz.eq.0) idz=dz
  if(iaa1.eq.-1) iaa1=aa1
  if(infuel.eq.0) infuel=1
endif
if (ifuelin.eq.0) then
  infuel=0
  inx=nx
  iny=ny
  inz=nz
  idx=dx
  idy=dy
  idz=dz
  iaa1=aa1
endif

!Find number species if using FastFuels dataset
if (itrees.eq.7) then
  call find_fastfuels_numspecies
  tfuelbins = 1
endif

end subroutine namelist_input

subroutine output_nfuel
!-----------------------------------------------------------------
! output is a function which writes the .dat files for use in 
! FIRETEC or QUIC-Fire
!-----------------------------------------------------------------
use grid_variables
use io_variables

implicit none

! Local Variables
integer ift,i,j,k,lfuel
real,dimension(nfuel):: nonzero
real*8,allocatable :: frhof(:,:,:,:),frhow(:,:,:,:)
real*8,allocatable :: fss(:,:,:,:),fafd(:,:,:)

! Executable Code
nonzero(:) = 1
lfuel = 1
do ift=1,nfuel
  if (sum(rhof(ift,:,:,:)).le.0) then
    nonzero(ift) = 0
  else
    do k=1,nz
      if (sum(rhof(ift,:,:,k)).gt.0) lfuel = max(lfuel,k)
    enddo
  endif
enddo

print*,'Exporting data to .dat files'
if (firetecshock.eq.1)then
  allocate(frhof(nfuel,nx,ny,lfuel),frhow(nfuel,nx,ny,lfuel), &
    fss(nfuel,nx,ny,lfuel),fafd(nfuel,nx,ny))
  frhof=rhof(:,:,:,1:lfuel)
  frhow=moist(:,:,:,1:lfuel)*rhof(:,:,:,1:lfuel)
  fss=sizescale(:,:,:,1:lfuel)
  fafd=fueldepth(:,:,:,1)
  open(1,file='fuelArrays.dat',form='unformatted',status='unknown')
  do ift=1,nfuel
    if (nonzero(ift).ne.0) then
      write (1) frhof(ift,:,:,:)
      write (1) frhow(ift,:,:,:)
      write (1) fss(ift,:,:,:)
      write (1) fafd(ift,:,:)
    endif
  enddo
  close (1)
  deallocate(frhof,frhow,fss,fafd)
else
  open (1,file='treesrhof.dat',form='unformatted',status='unknown')
  do ift=1,nfuel
    if (nonzero(ift).ne.0)  write (1) rhof(ift,:,:,:)
  enddo
  close (1)
  
  open (1,file='treesmoist.dat',form='unformatted',status='unknown')
  do ift=1,nfuel
    if (nonzero(ift).ne.0)  write (1) moist(ift,:,:,:)
  enddo
  close (1)
  
  open (1,file='treesss.dat',form='unformatted',status='unknown')
  do ift=1,nfuel
    if (nonzero(ift).ne.0)  write (1) sizescale(ift,:,:,:)
  enddo
  close (1)
  
  open (1,file='treesfueldepth.dat',form='unformatted',status='unknown')
  do ift=1,nfuel
    if (nonzero(ift).ne.0)  write (1) fueldepth(ift,:,:,:)
  enddo
  close (1)
endif

print*,'Your nfuel is',int(sum(nonzero(:)))
print*,'Your lfuel is',lfuel

end subroutine output_nfuel

subroutine output_1fuel
!-----------------------------------------------------------------
! output_1fuel is a function which writes the .dat files for use in 
! FIRETEC or QUIC-Fire but combines all fuels into a single fuel
! type.
!-----------------------------------------------------------------
use grid_variables
use io_variables

implicit none

! Local Variables
integer ift,i,j,k,lfuel
real,allocatable :: srhof(:,:,:),sss(:,:,:),smoist(:,:,:),safd(:,:,:)
real*8,allocatable :: frhof(:,:,:),frhow(:,:,:)
real*8,allocatable :: fss(:,:,:),fafd(:,:)

! Executable Code
allocate(srhof(nx,ny,nz))
allocate(sss(nx,ny,nz))
allocate(smoist(nx,ny,nz))
allocate(safd(nx,ny,nz))
do i=1,nx
  do j=1,ny
    do k=1,nz
      do ift=1,nfuel
        if(rhof(ift,i,j,k).gt.0)then
          sss(i,j,k)=(sss(i,j,k)*srhof(i,j,k)+sizescale(ift,i,j,k)*rhof(ift,i,j,k)) &
            /(srhof(i,j,k)+rhof(ift,i,j,k))
          smoist(i,j,k)=(smoist(i,j,k)*srhof(i,j,k)+moist(ift,i,j,k)*rhof(ift,i,j,k)) &
            /(srhof(i,j,k)+rhof(ift,i,j,k))
          safd(i,j,k)=(safd(i,j,k)*srhof(i,j,k)+fueldepth(ift,i,j,k)*rhof(ift,i,j,k)) &
            /(srhof(i,j,k)+rhof(ift,i,j,k))
          srhof(i,j,k)=srhof(i,j,k)+rhof(ift,i,j,k)
        endif
      enddo
    enddo
  enddo
enddo

lfuel = 1
do k=1,nz
  if (sum(srhof(:,:,k)).gt.0) lfuel = max(lfuel,k)
enddo

print*,'Exporting data to .dat files'
if (firetecshock.eq.1)then
  allocate(frhof(nx,ny,lfuel),frhow(nx,ny,lfuel), &
    fss(nx,ny,lfuel),fafd(nx,ny))
  frhof=srhof(:,:,1:lfuel)
  frhow=smoist(:,:,1:lfuel)*srhof(:,:,1:lfuel)
  fss=sss(:,:,1:lfuel)
  fafd=safd(:,:,1)
  open(1,file='fuelArrays.dat',form='unformatted',status='unknown')
  write (1) frhof(:,:,:)
  write (1) frhow(:,:,:)
  write (1) fss(:,:,:)
  write (1) fafd(:,:)
  close (1)
  deallocate(frhof,frhow,fss,fafd)
else
  open (1,file='treesrhof.dat',form='unformatted',status='unknown')
  write (1) srhof
  close (1)

  open (1,file='treesfueldepth.dat',form='unformatted',status='unknown')
  write (1) safd
  close (1)

  open (1,file='treesss.dat',form='unformatted',status='unknown')
  write (1) sss
  close (1)

  open (1,file='treesmoist.dat',form='unformatted',status='unknown')
  write (1) smoist
  close (1)
endif

print*,'Your nfuel is 1'
print*,'Your lfuel is',lfuel

end subroutine output_1fuel

subroutine find_fastfuels_numspecies
!-----------------------------------------------------------------
! This subroutine will find the number of species in the fast 
! fuels data
! Added 10/21 JO
!-----------------------------------------------------------------
use grid_variables
use baseline_variables
use species_variables

implicit none

!Local Variables
integer i,j,ift,ff_len,min_val_sp, max_val_sp
integer,allocatable :: final_uni_sp(:), uni_sp(:)
real,dimension(19):: temp_array ! FF trees csv has at least 19 columns

! Executable Code
ff_len = 0
open (2,file=treefile)
do
  read (2,*,end=19) !length of FF columns
  ff_len = ff_len+1
enddo
19  rewind(2)

allocate(tspecies(ff_len-1))
read(2,*) !read 1st line and throw away, has column headers
do i=1,ff_len-1
  read(2,*) temp_array(:)
  if (iFIAspecies.eq.1) then
    tspecies(i)=temp_array(1)
  else
    tspecies(i)=temp_array(5) !take from sp_grp, 5th pos
  !print*,'species# = ',temp_array(5)
  endif
enddo
rewind(2)

!find unique number of tree species!
allocate(uni_sp(maxval(tspecies)))
min_val_sp = minval(tspecies)-1
max_val_sp = maxval(tspecies)
i=0
do while (min_val_sp<max_val_sp)
  i = i+1
  min_val_sp = minval(tspecies, mask=tspecies>min_val_sp)
  uni_sp(i) = min_val_sp
enddo
allocate(final_uni_sp(i), source=uni_sp(1:i)) 
ntspecies = count(final_uni_sp==final_uni_sp)
print*,'New Set Species = ',ntspecies

deallocate(tspecies)
deallocate(uni_sp)
deallocate(final_uni_sp)

end subroutine find_fastfuels_numspecies


subroutine find_numspecies
!-----------------------------------------------------------------
! This subroutine will find the number of species in any treelist
! Added 10/22 JSM
!-----------------------------------------------------------------
use grid_variables
use baseline_variables
use species_variables

implicit none

!Local Variables
integer i,numtrees,numspec,min_val_sp, max_val_sp
integer,allocatable :: final_uni_sp(:), uni_sp(:)
integer,allocatable :: temp_array(:)

! Executable Code
numtrees = 0
open (2,file=treefile)
do
  read (2,*,end=5) !length of treelist columns
  numtrees = numtrees+1
enddo
5  rewind(2)

allocate(temp_array(numtrees))
do i=1,numtrees
   read(2,*) temp_array(i)
enddo

allocate(uni_sp(maxval(temp_array)))
min_val_sp = minval(temp_array)-1
max_val_sp = maxval(temp_array)

i=0

do while (min_val_sp<max_val_sp)
   i=i+1
   min_val_sp = minval(temp_array, mask=temp_array>min_val_sp)
   uni_sp(i) = min_val_sp
enddo
allocate(final_uni_sp(i), source=uni_sp(1:i))
numspec = count(final_uni_sp==final_uni_sp)

ntspecies = numspec

print*,'New Set Species Number = ',ntspecies

deallocate(temp_array)
deallocate(uni_sp)
deallocate(final_uni_sp)

end subroutine find_numspecies

