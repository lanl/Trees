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

subroutine writearrays

    use DUETio

    integer :: s,g

    !print*,'Subroutine writearrays running...'
    do s=1,domain%ns
        outarray%lrho(s,:,:) = sum(litter%lrho(s,:,:,:),DIM=3)
    enddo

    do g = 1,domain%ng
        outarray%frho(g,:,:,1) = outarray%srho(g,:,:) 
        outarray%fsss(g,:,:,1) = outarray%ssss(g,:,:) 
        outarray%fafd(g,:,:) = outarray%safd(g,:,:) 
        outarray%fh20(g,:,:,1) = outarray%sh20(g,:,:) 
    enddo
    
    do s = 1,domain%ns
        outarray%frho(domain%ng+domain%ns+s,:,:,:) = inarray%trhof(s,:,:,:)
        outarray%fh20(domain%ng+domain%ns+s,:,:,:) = inarray%moist(s,:,:,:)
        outarray%fafd(domain%ng+domain%ns+s,:,:) = inarray%depth(s,:,:)
        outarray%fsss(domain%ng+domain%ns+s,:,:,:) = inarray%sscale(s,:,:,:)
        
        outarray%srho(domain%ng+s,:,:) = outarray%lrho(s,:,:)
        outarray%ssss(domain%ng+s,:,:) = outarray%lsss(s,:,:)
        outarray%safd(domain%ng+s,:,:) = outarray%lafd(s,:,:)
        outarray%sh20(domain%ng+s,:,:) = outarray%lh20(s,:,:)

        outarray%frho(domain%ng+s,:,:,1) = outarray%srho(domain%ng+s,:,:) 
        outarray%fsss(domain%ng+s,:,:,1) = outarray%ssss(domain%ng+s,:,:)
        outarray%fafd(domain%ng+s,:,:) = outarray%safd(domain%ng+s,:,:) 
        outarray%fh20(domain%ng+s,:,:,1) = outarray%sh20(domain%ng+s,:,:) 
    enddo
    
    !print*,'Subroutine writearrays complete.'

end subroutine writearrays

!---------------------------------------------------------------------!

subroutine FF_writefiles

    use DUETio

    integer :: g,s,i,j

    real,allocatable :: surfrho(:,:,:),surfh20(:,:,:),surfafd(:,:,:),surfsss(:,:,:),trees(:,:)

    allocate(surfrho(domain%ns,domain%nx,domain%ny))
    allocate(surfh20(domain%ns,domain%nx,domain%ny))
    allocate(surfafd(domain%ns,domain%nx,domain%ny))
    allocate(surfsss(domain%ns,domain%nx,domain%ny))

    surfrho = outarray%srho
    surfh20 = outarray%sh20
    surfsss = outarray%ssss
    surfafd = outarray%safd

    !print*,'Subroutine FF_writefiles running...'

    open (12,file='surface_rhof_layered.dat',form='unformatted',status='unknown')
    write (12) surfrho !outarray%srho
    close (12)

    !print*,'1'

    open (12,file='surface_moist_layered.dat',form='unformatted',status='unknown')
    write (12) surfh20 !outarray%sh20
    close (12)

    !print*,'2'

    open (12,file='surface_ss_layered.dat',form='unformatted',status='unknown')
    write (12) surfsss !outarray%ssss
    close (12)

    !print*,'3'

    open (12,file='surface_depth_layered.dat',form='unformatted',status='unknown')
    write (12) surfafd !outarray%safd
    close (12)

    !print*,'4'

    open (12,file='surface_species.dat',form='formatted',status='unknown')
      do s=1,size(specarray)
        write (12,'(I16)') specarray(s)
      enddo
    close(12)

    !print*,'5'
  
    deallocate(surfrho)
    deallocate(surfh20)
    deallocate(surfafd)
    deallocate(surfsss)

    allocate(surfrho(2,domain%nx,domain%ny))
    !allocate(surfh20(2,domain%nx,domain%ny))
    !allocate(surfafd(2,domain%nx,domain%ny))
    !allocate(surfsss(2,domain%nx,domain%ny))

    !print*,'6'

    surfrho = 0.0
    !surfh20 = 0.0
    !surfafd = 0.0
    !surfsss = 0.0

    !print*,'7'

    do g=1,domain%ng
        surfrho(1,:,:) = surfrho(1,:,:) + outarray%srho(g,:,:)
        !surfh20(1,:,:) = surfh20(1,:,:) + outarray%sh20(g,:,:)
        !surfafd(1,:,:) = surfafd(1,:,:) + outarray%safd(g,:,:)
        !surfsss(1,:,:) = surfsss(1,:,:) + outarray%ssss(g,:,:)
    enddo

    !print*,'8'

    do s=1,domain%ns
        surfrho(2,:,:) = surfrho(2,:,:) + outarray%srho(domain%ng+s,:,:)
        !surfh20(2,:,:) = surfh20(2,:,:) + outarray%sh20(domain%ng+s,:,:)
        !surfafd(2,:,:) = surfafd(2,:,:) + outarray%safd(domain%ng+s,:,:)
        !surfsss(2,:,:) = surfsss(2,:,:) + outarray%ssss(domain%ng+s,:,:)
    enddo

    !print*,'9'

    open (12,file='surface_rhof.dat',form='unformatted',status='unknown')
    write (12) surfrho
    close (12)

    !open (12,file='surface_moist.dat',form='unformatted',status='unknown')
    !write (12) surfh20
    !close (12)
!
    !open (12,file='surface_ss.dat',form='unformatted',status='unknown')
    !write (12) surfsss
    !close (12)
!
    !open (12,file='surface_depth.dat',form='unformatted',status='unknown')
    !write (12) surfafd
    !close (12)

    !print*,'10'

    allocate(trees(domain%nx,domain%ny))

    trees = 0.0

    do i=1,domain%nx
        do j=1,domain%ny
            trees(i,j) = sum(inarray%trhof(:,i,j,:))
        enddo
    enddo


    open (12,file='flattrees.dat',form='unformatted',status='unknown')
    write (12) trees
    close (12)

    open (12,file='canopy.dat',form='unformatted',status='unknown')
    write (12) inarray%trhof
    close (12)

    !print*,'Subroutine FF_writefiles complete.'

    
end subroutine FF_writefiles

!---------------------------------------------------------------------!

subroutine TR_RunDUET(nx,ny,nz,ns,zheight,trhof,tmoist,tfueldepth,tsizescale,grhof,lrhof, &
    gmoist,lmoist,gsizescale,lsizescale,gfueldepth,lfueldepth)

    use DUETio
    use fillio
    use species
    use winds
    use mainLoop
    use support

    integer,intent(in) :: ns,nx,ny,nz

    real,intent(in) :: zheight(nx,ny,nz),tfueldepth(nx,ny,nz)
    real,intent(in) :: trhof(ns,nx,ny,nz),tmoist(ns,nx,ny,nz)
    real,intent(in) :: tsizescale(ns,nx,ny,nz)

    real,intent(inout) :: grhof(ns,nx,ny,nz),lrhof(ns,nx,ny,nz)
    real,intent(inout) :: gmoist(ns,nx,ny,nz),lmoist(ns,nx,ny,nz)
    real,intent(inout) :: gsizescale(ns,nx,ny,nz),lsizescale(ns,nx,ny,nz)
    real,intent(inout) :: gfueldepth(ns,nx,ny),lfueldepth(ns,nx,ny)

    integer :: g,l

    print*,'DUET running...'

    domain%nx = nx
    domain%ny = ny
    domain%nz = nz
    domain%ns = ns

    allocate(inarray%trhof(ns,nx,ny,nz))
    allocate(inarray%zheight(nx,ny,nz))
    allocate(inarray%moist(ns,nx,ny,nz))
    allocate(inarray%depth(ns,nx,ny))
    allocate(inarray%sscale(ns,nx,ny,nz))

    inarray%zheight(:,:,:)   = zheight(:,:,:)
    inarray%trhof(:,:,:,:)   = trhof(:,:,:,:)
    inarray%moist(:,:,:,:)   = tmoist(:,:,:,:) 
    inarray%depth(:,:,:)     = tfueldepth(:,:,:)
    inarray%sscale(:,:,:,:)  = tsizescale(:,:,:,:)

    call TR_filldomain
    call BuildSpeciesArray
    call makewinds
    call disperseLitter
    call decay
    call diffuse
    call depthMoistSs
    call growGrass
    call writearrays
    call FF_writefiles

    do g=1,domain%ng
        grhof(g,:,:,1) = outarray%srho(g,:,:)
        gmoist(g,:,:,1) = outarray%sh20(g,:,:)
        gsizescale(g,:,:,1) = outarray%ssss(g,:,:)
        gfueldepth(g,:,:) = outarray%safd(g,:,:)
    enddo

    do l=1,domain%ns
        lrhof(l,:,:,1) = outarray%lrho(l,:,:)
        lmoist(l,:,:,1) = outarray%lh20(l,:,:)
        lsizescale(l,:,:,1) = outarray%lsss(l,:,:)
        lfueldepth(l,:,:) = outarray%lafd(l,:,:)
    enddo
    print*,'Min and Max of litter density:'  ,minval(outarray%lrho),maxval(outarray%lrho)
    print*,'Min and Max of litter moisture:' ,minval(outarray%lh20),maxval(outarray%lh20)
    print*,'Min and Max of litter sizescale:',minval(outarray%lsss),maxval(outarray%lsss)
    print*,'Min and Max of litter depth:'    ,minval(outarray%lafd),maxval(outarray%lafd)
    print*,''
    print*,'Min and Max of grass density:'  ,minval(grhof),     maxval(grhof)    
    print*,'Min and Max of grass moisture:' ,minval(gmoist),    maxval(gmoist)  
    print*,'Min and Max of grass sizescale:',minval(gsizescale),maxval(gsizescale)
    print*,'Min and Max of grass depth:'    ,minval(gfueldepth),maxval(gfueldepth)

end subroutine TR_RunDUET




