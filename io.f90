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
   verbose

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

!if(iFIA.eq.1.and.itrees.ne.7) call find_numspecies
!allocate(FIA(ntspecies))

open(unit=15,file='fuellist',form='formatted',status='old')
     read (15,nml=fuellist)
     
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
     
     if (itrees.eq.0) then
       ntspecies=0
     endif
     
     if (igrass.eq.0) then
       duet_ngrass = ngrass
       ngrass=0
     endif

     if (itrees.eq.7) then
      if (iFIA.eq.1) then
        iFIAspecies=0
      end if
      call find_fastfuels_numspecies
      tfuelbins = 1
    endif
    
     if (infuel+ntspecies.eq.0) then
        iFIA = 0
        allocate(FIA(1*tfuelbins))        
     else
      allocate(FIA(infuel+ntspecies*tfuelbins))
      if (itrees.eq.7) then
        FIA = final_uni_sp
      else
        read (15,nml=speciesdata)
      end if
    endif
    
     if (ilitter.gt.0.and.ilitter.ne.2) then
       allocate(lrho(ntspecies*tfuelbins))
       allocate(lmoisture(ntspecies*tfuelbins))
       allocate(lss(ntspecies*tfuelbins))
       allocate(ldepth(ntspecies*tfuelbins))

       read (15,nml=litterdata)
     endif

     if (igrass.eq.1.or.ilitter.eq.2) then
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
tdnx(1)=floor(tdnx(1)/dx+1)
tdnx(2)=ceiling(tdnx(2)/dx)
tdny(1)=floor(tdny(1)/dy+1)
tdny(2)=ceiling(tdny(2)/dy)

!Find number species if using FastFuels dataset
if (itrees.eq.7) then
  
  
endif


if (iFIA.eq.1) then
  call define_species_variables
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
integer,allocatable :: uni_sp(:)
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
print*,'Groups = ',final_uni_sp

deallocate(tspecies)
deallocate(uni_sp)

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
integer,allocatable :: uni_sp(:)
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

end subroutine find_numspecies


!-------------------------------------------
!Write New Fuellist based on old params
!-------------------------------------------
subroutine output_fuellist
  use grid_variables
  use io_variables
  use infile_variables
  use baseline_variables
  use duet_variables
  use species_variables

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
  YearsSinceBurn,StepsPerYear

  namelist/speciesdata/ &
  FIA

  namelist/litterdata/ &
  lrho,lmoisture,lss,ldepth

  namelist/grassdata/ &
  grho,gmoisture,gss,gdepth,gmoistoverride

  namelist/winddata/ &
  uavg,vavg,ustd,vstd



  open(unit=15,file='fuellist',form='formatted',status='old')
    read (15,nml=fuellist)

    if (infuel+ntspecies.ne.0) then
       read (15,nml=speciesdata)
    end if

    if (ilitter.gt.0.and.ilitter.ne.2) then
      read (15,nml=litterdata)
    endif

    if (igrass.eq.1.or.ilitter.eq.2) then
      read (15,nml=grassdata)
    endif

    if (windprofile.eq.0) then
      read (15,nml=winddata)
    endif

  close(15)


  open(unit=2222,file='fuellist_'//newtreefile,form='formatted',status='unknown')
    !fuellist
    ndatax=nx
    ndatay=ny
    treefile=TRIM(newtreefile)//'_treelist.txt'
    write(2222,nml=fuellist)

    !species list
    if (infuel+ntspecies.ne.0) then
      write(2222,nml=speciesdata)
    end if

    if (ilitter.gt.0.and.ilitter.ne.2) then
      write (2222,nml=litterdata)
    endif

    if (igrass.eq.1.or.ilitter.eq.2) then
      write (2222,nml=grassdata)
    endif

    if (windprofile.eq.0) then
      read (2222,nml=winddata)
    end if


  close(2222)


end subroutine output_fuellist


!-------------------------------------------
!Write New Fuellist based on old params
!-------------------------------------------
subroutine output_1fuellist

use grid_variables
use io_variables
use infile_variables
use baseline_variables
use duet_variables
use species_variables

!stuff for writing new treelist
character(len=100) :: newtreelistname
integer :: ii

open (2222,file='fuellist_'//newtreefile,form='formatted',status='unknown')
newtreelistname = newtreefile//'_treelist.txt'
write(2222,'(A10)') '&fuellist'
write(2222,'(A35)') '! ----------------------------------'
write(2222,'(A35)') '! FIRETEC domain info'
write(2222,'(A35)') '! ----------------------------------'
write(2222,'(A4,I6.2,A5,I6.2,A5,I6.2,A50)') 'nx=',nx,' ny=',ny,' nz=',nz,'! Size of HIGRAD/FIRETEC grid [cells]'
write(2222,'(A4,F6.2,A5,F6.2,A5,F6.2,A50)') 'dx=',dx,' dy=',dy,' dz=',dz,'! Grid Resolution [m]'
write(2222,'(A4,F6.2,A50)') 'aa1=',aa1,'! Vertical stretching component [default=0.1]'
write(2222,'(A15,I2.1,A50)') 'singlefuel=',singlefuel,'! Flag forcing single fuel type instead of multiple fuels'
write(2222,'(3A30)') 'topofile= "'//TRIM(topofile)//'"  	     '
write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A35)')'! Data import from existing files'
write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A25)')'ifuelin=0 ! Fuel read-in flag '
write(2222,'(A50)')'inx=0 iny=0 inz=0 ! Size of HIGRAD/FIRETEC grid'
write(2222,'(A40)')'idx=0 idy=0 idz=0 ! Grid Resolution [m]'
write(2222,'(A50)')'iaa1=0.1 ! Vertical stretching component'
write(2222,'(A30)')'infuel=1 ! Number of Fuel Types'
write(2222,'(A15)')'intopofile= " " '
write(2222,'(A15)')'rhoffile =  " " '
write(2222,'(A15)')'ssfile   =  " " '
write(2222,'(A15)')'moistfile=  " " '
write(2222,'(A15)')'afdfile  =  " " '

write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A35)')'! Input dataset info'
write(2222,'(A35)')'! ----------------------------------'
if (itrees.eq.0) then 
  write(2222,'(A25)')'itrees=0 ! Trees flag '
else
  write(2222,'(A25)')'itrees=2 ! Trees flag '
end if
write(2222,'(A11,I6.3,A25)')'ntspecies=',ntspecies,'! Number of Tree Species'
write(2222,'(A5,I1,A25)')'iFIA=',iFIA,'!Turn on FIA code'
write(2222,'(A15,I1,A25)')'iFIAspecies=',iFIAspecies,'!species or groups '

write(2222,'(A12)') "tfuelbins=1"
write(2222,*) "treefile='",TRIM(newtreelistname)//'_treelist.txt',"'"
!write(2222,'(A10,A30,A1)') "treefile='",TRIM(newtreelistname),"'"
write(2222,'(A10,F6.1)') "ndatax=",nx*dx
write(2222,'(A10,F6.1)') "ndatay=",ny*dy
write(2222,'(A10)') "datalocx= "
write(2222,'(A10)') "datalocy= "
write(2222,'(A10)') "tdnx=     "
write(2222,'(A10)') "tdny=     "

write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A35)')'! Litter and grass settings'
write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A25,I1)')'ilitter=',ilitter
write(2222,'(A25,F4.2)')'litterconstant=',litterconstant
write(2222,'(A25,I1)')'windprofile=',windprofile
write(2222,'(A25,I2.1)')'randomwinds=',randomwinds
write(2222,'(A25,I3.1)')'grassstep=',grassstep
write(2222,'(A25,I4.1)')'YearsSinceBurn=',YearsSinceBurn
write(2222,'(A25,I4.1)')'StepsPerYear=',StepsPerYear
write(2222,'(A25,F6.2)')'relhum =',relhum
write(2222,'(A25,I1)')'igrass=',igrass
write(2222,'(A25,I2.1)')'ngrass=',ngrass
write(2222,'(A25,F4.2)')'grassconstant=',grassconstant

write(2222,'(A35)') '! ----------------------------------'
write(2222,'(A35)') '! Extra Options for '
write(2222,'(A35)') '! producing more output files and '
write(2222,'(A35)') '! controlling the random seed'
write(2222,'(A35)') '! ----------------------------------'
write(2222,'(A25,I4.1)')'litout =',litout
write(2222,'(A25,I4.1)')'controlseed =',controlseed
write(2222,'(A25,I4.1)')'seedchange =',seedchange

write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A50)')'! Option to output tree list and fuellist'
write(2222,'(A35)')'! ----------------------------------'
write(2222,'(A75)')'verbose = 0 !flag to output treelist and fuellist from resulting run: 0=No ; 1=Yes'
write(2222,'(A1)')'/'

!SPECIES NAMELIST=====================================================
write(2222,'(A15)')'&speciesdata'
write(2222,'(A5)', Advance = 'No' )'FIA ='
do ii = 1, ntspecies
  write(2222, '( I4.3, A1 )', Advance = 'No' ) FIA(ii)
  if (ii.ne.ntspecies) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/,A1)')'/'

!LITTER NAMELIST=====================================================
write(2222,'(A15)') '&litterdata'
write(2222,'(A10)', Advance = 'No' ) 'lrho ='
do ii = 1, ntspecies
  write(2222, '( F5.3, A1 )', Advance = 'No' ) lrho(ii)
  if (ii.ne.ntspecies) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/)')
write(2222,'(A15)', Advance = 'No' ) 'lmoisture =' 
do ii = 1, ntspecies
  write(2222, '( F4.3, A1 )', Advance = 'No' ) lmoisture(ii)
  if (ii.ne.ntspecies) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/)')
write(2222,'(A5)', Advance = 'No' ) 'lss =' 
do ii = 1, ntspecies
  write(2222, '( F6.5, A1 )', Advance = 'No' ) lss(ii)
  if (ii.ne.ntspecies) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/)')
write(2222,'(A10)', Advance = 'No' ) 'ldepth =' 
do ii = 1, ntspecies
  write(2222, '( F4.3, A1 )', Advance = 'No' ) ldepth(ii)
  if (ii.ne.ntspecies) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/,A1)')'/'


!GRASS NAMELIST=====================================================
write(2222,'(A15)') '&grassdata'

write(2222,'(A10)', Advance = 'No' ) 'grho ='
do ii = 1, ngrass
  write(2222, '( F5.4, A1 )', Advance = 'No' ) grho(ii) 
  if (ii.ne.ngrass) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/)')
write(2222,'(A11)', Advance = 'No' ) 'gmoisture ='
do ii = 1, ngrass
  write(2222, '( F4.3, A1 )', Advance = 'No' ) gmoisture(ii) 
  if (ii.ne.ngrass) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/)')
write(2222,'(A10)', Advance = 'No' ) 'gss ='
do ii = 1, ngrass
  write(2222, '( F6.4, A1 )', Advance = 'No' ) gss(ii)  
  if (ii.ne.ngrass) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/)')
write(2222,'(A10)', Advance = 'No' ) 'gdepth ='
do ii = 1, ngrass
  write(2222, '( F4.3, A1 )', Advance = 'No' ) gdepth(ii) 
  if (ii.ne.ngrass) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/)')
write(2222,'(A25,F3.1)') 'gmoistoverride = ',gmoistoverride
write(2222,'(A1)')'/'


!WIND NAMELIST=====================================================
if (windprofile.eq.0) then
  periodTotal=YearsSinceBurn*StepsPerYear
else
  periodTotal=0
end if 

write(2222,'(A15)') '&winddata'
write(2222,'(A10)', Advance = 'No' ) 'uavg = '
do ii = 1, periodTotal
  write(2222, '( F4.3, A1 )', Advance = 'No' ) uavg(ii) 
  if (ii.ne.periodTotal) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/)')
write(2222,'(A10)', Advance = 'No' ) 'vavg = '
do ii = 1, periodTotal
  write(2222, '( F4.3, A1 )', Advance = 'No' ) vavg(ii) 
  if (ii.ne.periodTotal) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/)')
write(2222,'(A10)', Advance = 'No' ) 'ustd = '
do ii = 1, periodTotal
  write(2222, '( F4.3, A1 )', Advance = 'No' ) ustd(ii)  
  if (ii.ne.periodTotal) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/)')
write(2222,'(A10)', Advance = 'No' ) 'vstd = '
do ii = 1, periodTotal
  write(2222, '( F4.3, A1 )', Advance = 'No' ) vstd(ii)  
  if (ii.ne.periodTotal) write(2222,'(A1)', Advance = 'No') ','
end do
write(2222,'(/,A1)')'/'

close(2222)
  

end subroutine output_1fuellist
