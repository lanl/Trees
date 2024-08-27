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

module FFio
    implicit none

    type :: FFinarrays
        real,allocatable :: FFrhof(:,:,:),FFmoist(:,:,:)
        integer*2,allocatable :: FFspec(:,:,:)
    end type FFinarrays

    type(FFinarrays) :: FFinarray

end module FFio

!---------------------------------------------------------------------!

module FastFuels

    contains

    subroutine FF_Readfiles

        use DUETio, only : domain,duetvars
        use FFio

        !integer :: seedchange,winddirection,windvary
        integer,allocatable :: seed(:)
        
        print*,'Subroutine ReadFFfiles running...'
        
        open(unit=1,file='duet.in',form='formatted',status='old')
        read(1,*) domain%nx
        read(1,*) domain%ny
        read(1,*) domain%nz
        read(1,*) domain%dx
        read(1,*) domain%dy
        read(1,*) domain%dz
        read(1,*) duetvars%seedchange
        read(1,*) duetvars%winddirection
        read(1,*) duetvars%windvary
        read(1,*) domain%YSB
        close(1)

        domain%SPY=1        
        domain%nt = domain%YSB*domain%SPY
        domain%ng = 1
        
        call random_seed(size=n)
        allocate(seed(n))
        seed = seedchange
        call random_seed(put=seed)
        deallocate(seed)
      
        allocate(FFinarray%FFrhof(domain%nx,domain%ny,domain%nz))
        allocate(FFinarray%FFmoist(domain%nx,domain%ny,domain%nz))
        allocate(FFinarray%FFspec(domain%nx,domain%ny,domain%nz))
        
        open(1,file='treesrhof.dat',form='unformatted',status='old')
        read(1) FFinarray%FFrhof
        close(1)
        
        open(2,file='treesmoist.dat',form='unformatted',status='old')
        read(2) FFinarray%FFmoist
        close(2)
        
        open(unit=1,file='treesspcd.dat',form='unformatted',status='unknown')
        read(1) FFinarray%FFSpec
        close(1)


        print*,'Subroutine ReadFFfiles complete.'
        
    end subroutine FF_Readfiles

    !---------------------------------------------------------------------!

    subroutine FF_Species

        use DUETio
        use FFio
        use fillio
    
        integer :: k,j,i,n
        integer*2,allocatable :: specs(:)
        !integer(2),intent(inout) :: FFinarray%FFSpec(domain%nx,domain%ny,domain%nz)
    
        integer,dimension(290) :: FIA
        integer,dimension(290) :: vals

        print*,'Subroutine FF_Species running...'
    
        FIA = 0
        vals = 0
    
        do i=1,domain%nx
          do j=1,domain%ny
            do k=1,domain%nz
              FIA(int(FFinarray%FFSpec(i,j,k))) = FIA(int(FFinarray%FFSpec(i,j,k))) + 1
            enddo
          enddo
        enddo
    
        n=1
        do i=1,290
          if (FIA(i).ne.0)  then 
            vals(n) = i
            n = n+1
          endif
        enddo
        
        allocate(specs(n-1))
        specs(:) = vals(1:n-1)
    
        domain%ns = int(n-1)

        call alloc_init
        allocate(inarray%trhof(domain%ns,domain%nx,domain%ny,domain%nz))
        allocate(inarray%moist(domain%ns,domain%nx,domain%ny,domain%nz))
        allocate(inarray%depth(domain%ns,domain%nx,domain%ny))
        allocate(inarray%sscale(domain%ns,domain%nx,domain%ny,domain%nz))
        allocate(inarray%zheight(domain%nx,domain%ny,domain%nz))

        specarray = specs

        print*,'specarray = ',specarray
        print*,'specs = ',specs

        grasses%decy  = 1.0
        grasses%fh20  = 0.06
        grasses%step  = 1
        grasses%grho  = 1.1766
        grasses%dept  = 0.27
        grasses%ssss  = 0.0005

        print*,'Subroutine FF_Species complete.'
    
        !print*,'specarray: ',specarray, vals, n
      
    end subroutine FF_Species

    !---------------------------------------------------------------------!
  
    subroutine FF_trhof

      use DUETio, only : domain,inarray,outarray,specarray
      use FFio
    
      integer :: s,k,j,i,z

      print*,'Subroutine FF_trhof running...'
      !print*,'domain',domain%ns,domain%nz,domain%ny,domain%nx
      !print*,'Shape of FFinarray%FFspec',shape(FFinarray%FFspec)
      !print*,'Shape of specarray',shape(specarray)
      !print*,'Shape of FFinarray%FFrhof',shape(FFinarray%FFrhof)
      !print*,'Shape of FFinarray%FFmoist',shape(FFinarray%FFmoist)
      !print*,'Shape of inarray%trhof',shape(inarray%trhof)
      !print*,'Shape of inarray%moist',shape(inarray%moist)
    
        do s = 1,domain%ns
          do k = 1,domain%nz
            do j = 1,domain%ny
              do i = 1,domain%nx
                if (FFinarray%FFspec(i,j,k).eq.specarray(s)) then

                  inarray%trhof(s,i,j,k) = FFinarray%FFrhof(i,j,k)
                  inarray%moist(s,i,j,k) = FFinarray%FFmoist(i,j,k)
    
                  outarray%frho(domain%ng+domain%ns+s,i,j,k) = inarray%trhof(s,i,j,k)
                  outarray%fh20(domain%ng+domain%ns+s,i,j,k) = inarray%moist(s,i,j,k)

                endif
              enddo
            enddo
          enddo
        enddo
        !print*,'WE GOT TO THIS POINT...'
        do z = 1,domain%nz
          inarray%zheight(:,:,z) = domain%dz*(z-1)
        enddo
        
        print*,'Max and Min of trhof:',maxval(inarray%trhof),minval(inarray%trhof)

        print*,'Subroutine FF_trhof complete.'
    
    end subroutine FF_trhof
    
end module FastFuels
