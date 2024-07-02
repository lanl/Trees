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
use baseline_variables
use duet_variables, only : windprofile, &
  grassstep,YearsSinceBurn,StepsPerYear,randomwinds,relhum, &
  litout,inputprogram,densitythresh
use species_variables
use FF_variables
implicit none

real :: PI

!Local Variables
namelist/duetlist/ &
   inputprogram, &
   nx,ny,nz,dx,dy,dz,zmax,PI, &
   fueltotal,infuel,singlefuel, &
   ntspecies,tfuelbins, &
   grassconstant,litterconstant, &
   ngrass,grassstep,gmoistoverride, &
   densitythresh,&
   iFIA,iFIAspecies, &
   windprofile,randomwinds, &
   winddirection,windvary, &
   YearsSinceBurn,StepsPerYear, &
   relhum,FuelFile,grassfile, &
   litout,controlseed,seedchange

 namelist/grassdata/ &
   grho,gmoisture,gss,gdepth

! Executable Code

open(unit=10,file='DUETinputs',form='formatted',status='old')
  read (10,nml=duetlist)
  
  allocate(grho(ngrass),gmoisture(ngrass),gss(ngrass),gdepth(ngrass))
  
  if(inputprogram.eq.1) then
    open(unit=101,file='grass_DMSR.dat',form='unformatted',status='unknown')
    read(101) gdepth,gmoisture,gss,grho
    close(101)
  elseif(inputprogram.eq.2) then
    read(10,nml=grassdata)
    print*,'Grass data = ',grho,gmoisture,gss,gdepth
  endif
close(10)

if (controlseed.eq.1) then
  call random_seed(size=n)
  allocate(seed(n))
  seed=seedchange
  call random_seed(put=seed)
  deallocate(seed)
endif

end subroutine namelist_input

subroutine output_nfuel
!-----------------------------------------------------------------
! output is a function which writes the .dat files for use in 
! FIRETEC or QUIC-Fire
!-----------------------------------------------------------------
use grid_variables
use io_variables
use baseline_variables

implicit none
 
! Local Variables
integer ift,k,lfuel
real,dimension(ngrass+2*fueltotal):: nonzero
real*8,allocatable :: frhof(:,:,:,:),frhow(:,:,:,:)
real*8,allocatable :: fss(:,:,:,:),fafd(:,:,:)

! Executable Code
nonzero(:) = 1
lfuel = 1
do ift=1,ngrass+2*fueltotal
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
  allocate(frhof(fueltotal,nx,ny,lfuel),frhow(fueltotal,nx,ny,lfuel), &
    fss(fueltotal,nx,ny,lfuel),fafd(fueltotal,nx,ny))
  frhof=rhof(:,:,:,1:lfuel)
  frhow=moist(:,:,:,1:lfuel)*rhof(:,:,:,1:lfuel)
  fss=sizescale(:,:,:,1:lfuel)
  fafd=fueldepth(:,:,:,1)
  open(1,file='fuelArrays.dat',form='unformatted',status='unknown')
  do ift=1,fueltotal
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
  do ift=1,2*fueltotal+ngrass
    if (nonzero(ift).ne.0)  write (1) rhof(ift,:,:,:)
  enddo
  close (1)
  print*,'Shape of rhof = ',shape(rhof)
  
  open (1,file='treesmoist.dat',form='unformatted',status='unknown')
  do ift=1,2*fueltotal+ngrass
    if (nonzero(ift).ne.0)  write (1) moist(ift,:,:,:)
  enddo
  close (1)
  
  open (1,file='treesss.dat',form='unformatted',status='unknown')
  do ift=1,2*fueltotal+ngrass
    if (nonzero(ift).ne.0)  write (1) sizescale(ift,:,:,:)
  enddo
  close (1)
  
  open (1,file='treesfueldepth.dat',form='unformatted',status='unknown')
  do ift=1,2*fueltotal+ngrass
    if (nonzero(ift).ne.0)  write (1) fueldepth(ift,:,:,:)
  enddo
  close (1)
endif

print*,'Your fueltotal is',int(sum(nonzero(:)))
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
use baseline_variables

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
      do ift=1,fueltotal
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

print*,'Your fueltotal is 1'
print*,'Your lfuel is',lfuel

end subroutine output_1fuel


subroutine output_FF

use FF_variables, only : surfrhof,surfdepth,specarray

implicit none

integer :: s

! Executable Code

open (11,file='surface_rhof.dat',form='unformatted',status='unknown')
write (11) surfrhof
close (11)

open (111,file='surface_depth.dat',form='unformatted',status='unknown')
write (111) surfdepth
close (111)

open (1111,file='surface_species.dat',form='formatted',status='unknown')
!write (1111) specarray
  do s=1,size(specarray)
    print*,'Verifying specarray:',specarray(s)
    write (1111,'(I16)') specarray(s)
  enddo
close(1111)


end subroutine output_FF

