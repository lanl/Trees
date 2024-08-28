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


module species
contains

  subroutine BuildSpeciesArray

    use DUETio

    integer :: s,f

    !print*,'Subroutine BuildSpeciesArray running...'

    if(duetvars%iFIA.eq.1) then
      if(duetvars%iFIAspecies.eq.1) then
        do s=1,domain%ns
          do f=1,290
            if(SPECINFO(f)%FIA_code.eq.specarray(s)) then
              species%fuSA(s)=SPECINFO(f)%surfarea
              species%drop(s)=SPECINFO(f)%dropperyear
              species%decy(s)=SPECINFO(f)%decay
              species%step(s)=int(ceiling(real(SPECINFO(f)%stepperyear)/(12/domain%spy)))
              species%vter(s)=SPECINFO(f)%vterminal
              species%frde(s)=SPECINFO(f)%froude
              species%drag(s)=SPECINFO(f)%dragco
              species%fh20(s)=SPECINFO(f)%moist
              species%ssss(s)=SPECINFO(f)%sizescale
              species%comp(s)=SPECINFO(f)%compact
              exit
            endif
          enddo
        enddo
      else 
        do s=1,domain%ns
          do f=1,10
            if(f.eq.specarray(s)) then
              species%fuSA(s)=SPECgroups(f)%surfarea
              species%drop(s)=SPECgroups(f)%dropperyear
              species%decy(s)=SPECgroups(f)%decay
              species%step(s)=int(ceiling(real(SPECgroups(f)%stepperyear)/(12/domain%spy)))
              species%vter(s)=SPECgroups(f)%vterminal
              species%frde(s)=SPECgroups(f)%froude
              species%drag(s)=SPECgroups(f)%dragco
              species%fh20(s)=SPECgroups(f)%moist
              species%ssss(s)=SPECgroups(f)%sizescale
              species%comp(s)=SPECgroups(f)%compact
              exit
            endif
          enddo
        enddo
      endif
    endif

    !print*,'Subroutine BuildSpeciesArray complete.'

  end subroutine BuildSpeciesArray

  !---------------------------------------------------------------------!

  subroutine TR_species
    use DUETio

    integer :: s

    integer,dimension(4) :: trhofshape

    !print*,'Subroutine TR_species running...'

    trhofshape = shape(inarray%trhof)

    domain%ns = trhofshape(1)

    do s=1,domain%ns
      outarray%frho(domain%ng+domain%ns+s,:,:,:) = inarray%trhof(s,:,:,:)
    enddo

    !print*,'Subroutine TR_species complete.'
  end subroutine TR_species

end module species

!---------------------------------------------------------------------!

module winds
  implicit none
  contains
  subroutine makewinds

    use DUETio, only : windarray,duetvars,domain

    integer :: yt
    !integer :: windprofile,randomwinds
    real :: low,high,a
    real,allocatable :: Umean(:,:,:),Vmean(:,:,:),Uvar(:,:,:),Vvar(:,:,:)

    !print*,'Subroutine makewinds running...'

    if(duetvars%windprofile.eq.0) then
    
      print*,'!----!----!----!----!----!----!----!----!----!----!----!'
      print*,'Windprofile is taken from a user provided in fuellist'
      open(3,file='windprofile.dat',form='unformatted',status='old')
        read(3) windarray%uavg,windarray%vavg,windarray%uvar,windarray%vvar
      close(3)
    
      elseif(duetvars%windprofile.eq.1)then
      print*,'!----!----!----!----!----!----!----!----!----!----!----!'
      print*,'Windprofile is generated from higrad and averaged within each cell'
      
      allocate(Umean(domain%nx,domain%ny,domain%nz))
      allocate(Vmean(domain%nx,domain%ny,domain%nz))
      allocate(Uvar(domain%nx,domain%ny,domain%nz))
      allocate(VVar(domain%nx,domain%ny,domain%nz))
      
      open (105,file='Umean.dat', form='unformatted', status='old')
      read (105) Umean
      close(105)
    
      open (106,file='Vmean.dat', form='unformatted', status='old')
      read (106) Vmean
      close(106)
    
      open (107,file='Uvar.dat', form='unformatted', status='old')
      read (107) Uvar
      close(107)
    
      open (108,file='Vvar.dat', form='unformatted', status='old')
      read (108) Vvar
      close(108)
    
    elseif(duetvars%windprofile.eq.2) then
      print*,'!----!----!----!----!----!----!----!----!----!----!----!' 
      print*,'Windprofile is randomly generated averages over whole domain' 
      print*,'multiplication factor for winds = ',duetvars%randomwinds,'over',domain%nt,'time periods'
      low=1*duetvars%randomwinds
      high=3*duetvars%randomwinds
      
      do yt=1,domain%nt
        a = rand()
        a = 2.0*a - 1.0
        if (a.gt.0) then
          windarray%uavg(yt) = (a*real(low+high))-real(low)
        else 
          windarray%uavg(yt) = (a*real(low+high))+real(low)
        endif
        a = rand()
        a = 2.0*a - 1.0
        if (a.gt.0) then
          windarray%vavg(yt) = (a*real(low+high))-real(low)
        else
          windarray%vavg(yt) = (a*real(low+high))+real(low)
        endif
        a = rand()    
        a = 2.0*a - 1.0
        if (a.gt.0) then
          windarray%uvar = (a*real(low+high))-real(low)
        else
          windarray%uvar = (a*real(low+high))+real(low)
        endif
        a = rand()     
        a = 2.0*a - 1.0
        if (a.gt.0) then
          windarray%vvar = (a*real(low+high))-real(low)
        else
          windarray%vvar = (a*real(low+high))+real(low)
        endif
      enddo
    
    endif
    print*,'uavg = ',windarray%uavg
    print*,'vavg = ',windarray%vavg
    print*,'uvar = ',windarray%uvar
    print*,'vvar = ',windarray%vvar
    print*,'Each column is a year'
    print*,'!----!----!----!----!----!----!----!----!----!----!----!' 
    
    !print*,'Subroutine makewinds complete.'

  end subroutine makewinds

end module winds
    