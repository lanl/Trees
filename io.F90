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

!-----------------------------------------------------------------------
! namelist_input is a function reads in the user-defined variables 
! from the treelist then assigns those variables for use by the 
! trees softeware.
!-----------------------------------------------------------------------
subroutine fuellist_input
  use grid_variables, only : nx,ny,nz,dx,dy,dz,aa1,topofile,ndatax, &
    ndatay,datalocx,datalocy
  use io_variables, only : workdir,filesep,singlefuel,lreduced, &
    doubleprec,verbose
  use infile_variables, only : ifuelin,inx,iny,inz,idx,idy,idz,iaa1, &
    infuel,intopofile,rhoffile,ssfile,moistfile,afdfile
  use fuel_variables, only : ilitter,igrass,ngrass,grassconstant, &
    grho,gmoisture,gss,gdepth,itrees,tfuelbins,treefile,ntspecies, &
    trhofmax,litterconstant,lrho,lmoisture,lss,ldepth
  implicit none
  
  ! Local Variables
  
  ! Executable Code
  !Set working directory!
  call SetWorkingDirectory(workdir,filesep)
  print*,workdir
  print*,filesep
  print*,TRIM(TRIM(workdir)//filesep)//'fuellist'
  open(unit=48,file=TRIM(TRIM(workdir)//filesep)//'fuellist', &
    form= 'formatted',status= 'old')
  
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
    call QueryFuellist_string('rhoffile',rhoffile,48, &
      'treesrhof.dat.orig')
    call QueryFuellist_string('ssfile',ssfile,48,'treesss.dat.orig')
    call QueryFuellist_string('moistfile',moistfile,48, &
      'treesmoist.dat.orig')
    call QueryFuellist_string('afdfile',afdfile,48, &
      'treesfueldepth.dat.orig')
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
    call QueryFuellist_string('treefile',treefile,48, &
      'AJoseEglinTrees.txt')
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
    call QueryFuellist_real_array('lrho',lrho,ntspecies*tfuelbins, &
      48,1.18)
    allocate(lmoisture(ntspecies*tfuelbins))
    call QueryFuellist_real_array('lmoisture',lmoisture, &
      ntspecies*tfuelbins,48,0.06)
    allocate(lss(ntspecies*tfuelbins))
    call QueryFuellist_real_array('lss',lss,ntspecies*tfuelbins, &
      48,0.0005)
    allocate(ldepth(ntspecies*tfuelbins))
    call QueryFuellist_real_array('ldepth',ldepth,ntspecies*tfuelbins, &
      48,0.05)
  endif
  
  ! Verbose
  call QueryFuellist_integer('verbose',verbose,48,0)
       
  close(48)

end subroutine fuellist_input

!-----------------------------------------------------------------------
! For use to specifiy a different working directory
! than the one where the exe is located
!-----------------------------------------------------------------------
subroutine SetWorkingDirectory(workingVariable,fileVariable)
  Implicit None
  
  ! Local Variables
  character(len=255),intent(inout) :: workingVariable
  character(len=1),intent(inout)   :: fileVariable
  
  integer :: i,narg,k
  character(255) :: instring
  
  ! Executable Code
  k = COMMAND_ARGUMENT_COUNT()
  print *,k
  IF(k > 0)THEN
    narg = 1
  ELSE
    narg = 0
  ENDIF
  
  CALL GETARG(narg, instring)
  if(narg == 1) THEN
    workingVariable = TRIM(instring)
    DO i = 1,256
      IF(instring(i:i) .eq. '/')THEN
        fileVariable = '/'
        EXIT
      ELSEIF(instring(i:i) .eq. '\')THEN
        fileVariable = '\'
        EXIT
      ENDIF
    ENDDO
    print *,'Working directory is set to: '//TRIM(workingVariable)
  else
    workingVariable = ''
    fileVariable = ''
  endif
end subroutine SetWorkingDirectory

!-----------------------------------------------------------------------
! output is a function which writes the .dat files for use in 
! FIRETEC or QUIC-Fire
!-----------------------------------------------------------------------
subroutine output_fuel
  use grid_variables, only : nx,ny,nz,rhof,moist,fueldepth,sizescale, &
    nfuel
  use io_variables, only : singlefuel,doubleprec,lreduced,workdir, &
    filesep
  implicit none
  
  ! Local Variables
  integer ift,i,j,k,lfuel,lz
  integer,dimension(nfuel):: nonzero
  real(8),allocatable :: rho(:,:,:,:),h2o(:,:,:,:)
  real(8),allocatable :: sss(:,:,:,:),afd(:,:,:,:)
  
  ! Executable Code
  if (singlefuel.eq.1) then ! Condense fuels to a single fuel type
    do i=1,nx
      do j=1,ny
        do k=1,nz
          if(sum(rhof(:,i,j,k)).gt.0)then
            sizescale(1,i,j,k)=sum(sizescale(:,i,j,k)*rhof(:,i,j,k))/ &
              sum(rhof(:,i,j,k))
            moist(1,i,j,k)=sum(moist(:,i,j,k)*rhof(:,i,j,k))/ &
              sum(rhof(:,i,j,k))
            fueldepth(1,i,j,k)=sum(fueldepth(:,i,j,k)*rhof(:,i,j,k))/ &
              sum(rhof(:,i,j,k))
            rhof(1,i,j,k)=sum(rhof(:,i,j,k))
            if(nfuel.gt.1) rhof(2:nfuel,i,j,k)=0.
          endif
        enddo
      enddo
    enddo
  endif
    
  ! Eliminate empty fuels
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
  if (lreduced.eq.1) then
    lz=lfuel
  else
    lz=nz
  endif
 
  ! Write the fuel files 
  print*,'Exporting data to .dat files'
  if(doubleprec.eq.1) then
    print*,'Converting to Double Precision'
    allocate(rho(nfuel,nx,ny,nz)); rho=dble(rhof)
    allocate(h2o(nfuel,nx,ny,nz)); h2o=dble(moist)
    allocate(sss(nfuel,nx,ny,nz)); sss=dble(sizescale)
    allocate(afd(nfuel,nx,ny,nz)); afd=dble(fueldepth)
  
    open (1,file=TRIM(TRIM(workdir)//filesep)//'treesrhof.dat', &
      form='unformatted',status='unknown')
    do ift=1,nfuel
      if (nonzero(ift).ne.0)  write (1) rho(ift,:,:,1:lz)
    enddo
    close (1)
    open (2,file=TRIM(TRIM(workdir)//filesep)//'treesmoist.dat', &
      form='unformatted',status='unknown')
    do ift=1,nfuel
      if (nonzero(ift).ne.0)  write (1) h2o(ift,:,:,1:lz)
    enddo
    close (2)
    open (3,file=TRIM(TRIM(workdir)//filesep)//'treesss.dat', &
      form='unformatted',status='unknown')
    do ift=1,nfuel
      if (nonzero(ift).ne.0)  write (1) sss(ift,:,:,1:lz)
    enddo
    close (3)
    open (4,file=TRIM(TRIM(workdir)//filesep)//'treesfueldepth.dat', &
      form='unformatted',status='unknown')
    do ift=1,nfuel
      if (nonzero(ift).ne.0)  write (1) afd(ift,:,:,1:lz)
    enddo
    close (4)

    deallocate(rho,h2o,sss,afd)
  else
    open (1,file=TRIM(TRIM(workdir)//filesep)//'treesrhof.dat', &
      form='unformatted',status='unknown')
    do ift=1,nfuel
      if (nonzero(ift).ne.0)  write (1) rhof(ift,:,:,1:lz)
    enddo
    close (1)
    open (2,file=TRIM(TRIM(workdir)//filesep)//'treesmoist.dat', &
      form='unformatted',status='unknown')
    do ift=1,nfuel
      if (nonzero(ift).ne.0)  write (2) moist(ift,:,:,1:lz)
    enddo
    close (2)
    open (3,file=TRIM(TRIM(workdir)//filesep)//'treesss.dat', &
      form='unformatted',status='unknown')
    do ift=1,nfuel
      if (nonzero(ift).ne.0)  write (3) sizescale(ift,:,:,1:lz)
    enddo
    close (3)
    open (4,file=TRIM(TRIM(workdir)//filesep)//'treesfueldepth.dat', &
      form='unformatted',status='unknown')
    do ift=1,nfuel
      if (nonzero(ift).ne.0)  write (4) fueldepth(ift,:,:,1:lz)
    enddo
    close (4)
  endif
  
  print*,'Min and Max of rhof layer 1:', &
    minval(rhof(:,:,:,1)),maxval(rhof(:,:,:,1))
  print*,'Your nfuel is',int(sum(nonzero(:)))
  print*,'Your lfuel is',lfuel
end subroutine output_fuel
  
!-----------------------------------------------------------------------
! Write New Fuellist based on old params
!-----------------------------------------------------------------------
subroutine output_fuellist
  use grid_variables, only : nx,ny,nz,dx,dy,dz,aa1,topofile, &
    datalocx,datalocy
  use infile_variables, only : ifuelin
  use io_variables, only : workdir,filesep,singlefuel
  use fuel_variables, only : ilitter,igrass,ngrass,grassconstant, &
    grho,gmoisture,gss,gdepth,itrees,tfuelbins,ntspecies, &
    litterconstant,lrho,lmoisture,lss,ldepth
  use fuel_variables, only : newtreefile
  Implicit None
  
  ! Local Variables
  character(len=300) :: newtreelistname
  integer :: ii
  
  open (2222,file= TRIM(TRIM(workdir)//filesep)//'fuellist_'// &
    newtreefile,form= 'formatted',status= 'unknown')
  newtreelistname = newtreefile//'_treelist.txt'
  write(2222,'(A35)') '! ----------------------------------'
  write(2222,'(A35)') '! FIRETEC domain info'
  write(2222,'(A35)') '! ----------------------------------'
  write(2222,'(A4,I6.2)')     'nx = ',nx
  write(2222,'(A4,I6.2)')     'ny = ',ny
  write(2222,'(A4,I6.2,A50)') 'nz = ',nz
  write(2222,'(A4,F6.2)')     'dx = ',dx
  write(2222,'(A4,F6.2)')     'dy = ',dy
  write(2222,'(A4,F6.2,A50)') 'dz = ',dz
  write(2222,'(A6,F6.2,A50)') 'aa1 = ',aa1
  write(2222,'(A15,I2.1,A60)') 'singlefuel = ',singlefuel
  write(2222,'(3A30)') "topofile = '"// &
    TRIM(workdir//filesep//topofile)//"' "
  
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
      write(2222,*) "treefile = '", &
        TRIM(workdir//filesep//newtreelistname)//'_treelist.txt',"'"
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
      if (ii.ne.ntspecies*tfuelbins) write(2222,'(A1)',Advance='No') ','
    enddo
    write(2222,'(/,A15)', Advance = 'No' ) 'lmoisture = ' 
    do ii = 1, ntspecies*tfuelbins
      write(2222, '( F4.3, A1 )', Advance = 'No' ) lmoisture(ii)
      if (ii.ne.ntspecies*tfuelbins) write(2222,'(A1)',Advance='No') ','
    enddo
    write(2222,'(/,A6)', Advance = 'No' ) 'lss = ' 
    do ii = 1, ntspecies*tfuelbins
      write(2222, '( F6.5, A1 )', Advance = 'No' ) lss(ii)
      if (ii.ne.ntspecies*tfuelbins) write(2222,'(A1)',Advance='No') ','
    enddo
    write(2222,'(/,A10)', Advance = 'No' ) 'ldepth = ' 
    do ii = 1, ntspecies*tfuelbins
      write(2222, '( F4.3, A1 )', Advance = 'No' ) ldepth(ii)
      if (ii.ne.ntspecies*tfuelbins) write(2222,'(A1)',Advance='No') ','
    enddo
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
  write(2222,'(A75)')'verbose = 0'
  
  close(2222)
    
end subroutine output_fuellist
