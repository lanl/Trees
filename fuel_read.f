      subroutine grid_readin
      !-----------------------------------------------------------------
      ! grid_readin is a function which reads in .dat files for use in 
      ! FIRETEC or QUIC-Fire
      !-----------------------------------------------------------------
      use grid_variables
      use infile_variables
      use baseline_variables
      implicit none
      integer ift,i,j,k
      integer ii,jj,kk,kkk
      integer xbot,xtop,ybot,ytop,zbot,ztop
      real cells,xfrac,yfrac,zfrac
      real zstemp,rhoftemp 
      real target_mass,actual_mass
      real,dimension(2):: xcor,ycor,zcor

      allocate(irhof(infuel,inx,iny,inz))
      allocate(imoist(infuel,inx,iny,inz))
      allocate(iss(infuel,inx,iny,inz))
      allocate(iafd(infuel,inx,iny,inz))

      ! Read in fuel files
      print*,"Reading existing fuel files"
      open(unit=1,file=rhoffile,form='unformatted',status='unknown')
        read(1) irhof
      close(1)
      where (irhof<0)
        irhof = 0
      endwhere
      print*,'irhof','min=',minval(irhof),'max=',maxval(irhof),'avg=',sum(irhof)/(inx*iny*inz)
      open(unit=2,file=moistfile,form='unformatted',status='unknown')
        read(2) imoist
      close(2)
      where (imoist<0)
        imoist = 0
      endwhere
      print*,'imoist','min=',minval(imoist),'max=',maxval(imoist),'avg=',sum(imoist)/(inx*iny*inz)
      open(unit=3,file=ssfile,form='unformatted',status='unknown')
        read(3) iss
      close(3)
      if (ssfile.eq."sav.dat.orig")then ! Special case if ss file was sav
        do i=1,inx
          do j=1,iny
            do k=1,inz
              do ift=1,infuel
                if(iss(ift,i,j,k).gt.1) iss(ift,i,j,k)=2./iss(ift,i,j,k)
              enddo
            enddo
          enddo
        enddo
      endif
      where (iss<0)
        iss = 0
      endwhere
      print*,'iss','min=',minval(iss),'max=',maxval(iss),'avg=',sum(iss)/(inx*iny*inz)
      open(unit=4,file=afdfile,form='unformatted',status='unknown')
        read(4) iafd
      close(4)
      where (iafd<0)
        iafd = 0
      endwhere
      print*,'iafd','min=',minval(iafd),'max=',maxval(iafd),'avg=',sum(iafd)/(inx*iny*inz)

      print*,"Exisiting fuel files readin"

      ! Interpolate read file onto FIRETEC grid
      if (iintpr.eq.1) then
        print*,"Interpolating readin fuel files to desired FIRETEC grid"
        do i=1,nx
          do j=1,ny
            xcor(1) = (i-1)*dx ! Real x lower edge
            xcor(2) = i*dx     ! Real x upper edge
            xbot    = floor(xcor(1)/idx+1) ! Fuel readin grid x lower edge
            xtop    = min(inx,floor(xcor(2)/idx+1)) ! Fuel readin grid x upper edge
            ycor(1) = (j-1)*dy ! Real y lower edge
            ycor(2) = j*dy     ! Real y upper edge
            ybot    = floor(ycor(1)/idy+1) ! Fuel readin grid y lower edge
            ytop    = min(iny,floor(ycor(2)/idy+1)) ! Fuel readin grid y upper edge
            do k=1,nz-1
              cells = 0.
              do ii=xbot,xtop
                do jj=ybot,ytop
                  do kk=1,inz-1
                    if(izheight(ii,jj,kk).ge.zheight(i,j,k+1)) exit
                    if(izheight(ii,jj,kk).ge.zheight(i,j,k))then
                      do ift=1,infuel
                        xfrac = (min(ii*idx,xcor(2))-max((ii-1)*idx,xcor(1)))/idx
                        yfrac = (min(jj*idy,ycor(2))-max((jj-1)*idy,ycor(1)))/idy
                        zfrac = (min(izheight(ii,jj,kk+1),zheight(i,j,k+1))-max(izheight(ii,jj,kk),zheight(i,j,k)))/
     +                    (izheight(ii,jj,kk+1)-izheight(ii,jj,kk))
                        cells = cells+xfrac*yfrac*zfrac
                        rhoftemp=irhof(ift,ii,jj,kk)*xfrac*yfrac*zfrac
                        if(rhoftemp.gt.0)then
                        sizescale(ift,i,j,k) = (rhof(ift,i,j,k)*sizescale(ift,i,j,k)+rhoftemp*iss(ift,ii,jj,kk))/
     +                    (rhof(ift,i,j,k)+rhoftemp)
                        moist(ift,i,j,k) = (rhof(ift,i,j,k)*moist(ift,i,j,k)+rhoftemp*imoist(ift,ii,jj,kk))/
     +                    (rhof(ift,i,j,k)+rhoftemp)
                        fueldepth(ift,i,j,k) = (rhof(ift,i,j,k)*fueldepth(ift,i,j,k)+rhoftemp*iafd(ift,ii,jj,kk))/
     +                    (rhof(ift,i,j,k)+rhoftemp)
                        rhof(ift,i,j,k)=rhof(ift,i,j,k)+rhoftemp
                        endif
                      enddo
                    endif
                  enddo
                enddo
              enddo
              do ift=1,infuel
                if(cells.gt.0) rhof(ift,i,j,k)=rhof(ift,i,j,k)/cells
              enddo
            enddo
          enddo
        enddo
      endif
      
      ! Print out the target and actual fuel masses for comparisons sake
      target_mass = 0
      do ift=1,infuel
        do i=1,min(inx,nint(nx/dx*idx))
          do j=1,min(iny,nint(ny/dy*idy))
            do k=1,inz-1
              target_mass = target_mass+irhof(ift,i,j,k)*idx*idy*(izheight(i,j,k+1)-izheight(i,j,k))
            enddo
          enddo
        enddo
      enddo
      print*,'Readin target fuel mass:',target_mass
      actual_mass = 0
      do ift=1,infuel
        do i=1,nx
          do j=1,ny
            do k=1,nz-1
              actual_mass = actual_mass+rhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
            enddo
          enddo
        enddo
      enddo
      print*,'Readin actual fuel mass:',actual_mass
      print*,'Readin error:',actual_mass/target_mass*100,'%'
       
      deallocate(irhof)
      deallocate(imoist)
      deallocate(iss)
      deallocate(iafd)
      
      end subroutine grid_readin
