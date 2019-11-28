      subroutine output_hg(plotvar,input)

! You should not need to touch this subroutine.

! This subroutine writes a file to unit 7 for use by the plotting
! program "eulplt".

      use common_block


      implicit none

! Local stuff
      integer ::   input, i,j

      real, dimension(i_max,j_max)  ::   pstag, vmach, plotvar
      real       eke,tstat,velsq,tstag

      open(unit=7,file='euler.plt')

      write(7,100) title
100   format(a80)

      cp = rgas*gamma/(gamma-1.)
      cv = cp/gamma
      write(7,101) 1,ni,nj,0,0,1,ni,0,0
      write(7,102) cp,gamma

      do i=1,ni
        write(7,102) (x(i,j),j=1,nj)
        write(7,102) (y(i,j),j=1,nj)
      end do
101   format(16i5)
102   format(10f10.5)

! Calculate the secondary variables

      if( input==1 ) then
        do i=1,ni
          do j=1,nj
            vx(i,j)  = rovx(i,j)/ro(i,j)
            vy(i,j)  = rovy(i,j)/ro(i,j)
            eke      = 0.5*(vx(i,j)*vx(i,j) + vy(i,j)*vy(i,j))
            tstat    = (hstag(i,j) - eke)/cp
            p(i,j)   = ro(i,j)*rgas*tstat
            roe(i,j) = ro(i,j)*(cv*tstat + eke)
          end do
        end do

! Calculate the mach number and stagnation presssure

        do i=1,ni
          do j=1,nj
            tstat = p(i,j)/rgas/ro(i,j)
            velsq = vx(i,j)*vx(i,j) + vy(i,j)*vy(i,j)
            tstag = tstat + 0.5*velsq/cp
            vmach(i,j) = sqrt(velsq/(gamma*rgas*tstat))
            pstag(i,j) = p(i,j)*(tstag/tstat)**(1/fga)
          end do
        end do

        write(7,*) ' time step number', nstep

        write(7,*)  ' axial velocity'
        do i=1,ni
          write(7,103) (vx(i,j),j=1,nj)
        end do
103     format(10f10.4)

        write(7,*)  ' y  velocity'
        do i=1,ni
          write(7,103) (vy(i,j),j=1,nj)
        end do

        write(7,*)  ' radial velocity'
        do i=1,ni
          write(7,103) (0.0,j=1,nj)
        end do

        write(7,*)  ' mach number '
        do i=1,ni
          write(7,103) (vmach(i,j),j=1,nj)
        end do

        write(7,*)  ' static pressure'
        do i=1,ni
          write(7,104) (p(i,j),j=1,nj)
        end do
104     format(10f10.1)

        write(7,*)  ' density '
        do i=1,ni
          write(7,103) (ro(i,j),j=1,nj)
        end do

        write(7,*)  ' variable plotvar '
        do i=1,ni
          write(7,103) (plotvar(i,j),j=1,nj)
        end do

        write(7,*)  ' del rovx '
        do i=1,ni
          write(7,105) (delrovx(i,j),j=1,nj)
        end do
105     format(10e10.4)

! This is the dummy part used only for initial stages of grid debugging
      else

        nstep = -1
        do i=1,ni
          do j=1,nj
            ro(i,j)  = 1.226
            vx(i,j)  = 0.
            vy(i,j)  = 0.
            tstat    = 288.15
            p(i,j)   = ro(i,j)*rgas*tstat
            roe(i,j) = ro(i,j)*cv*tstat
          end do
        end do

! Calculate the mach number and stagnation presssure

        do i=1,ni
          do j=1,nj
            tstat = p(i,j)/rgas/ro(i,j)
            velsq = vx(i,j)*vx(i,j) + vy(i,j)*vy(i,j)
            tstag = tstat + 0.5*velsq/cp
            vmach(i,j) = sqrt(velsq/(gamma*rgas*tstat))
            pstag(i,j) = p(i,j)*(tstag/tstat)**(1/fga)
          end do
        end do

        write(7,*) ' time step number', nstep

        write(7,*)  ' axial velocity'
        do i=1,ni
          write(7,103) (vx(i,j),j=1,nj)
        end do

        write(7,*)  ' y  velocity'
        do i=1,ni
          write(7,103) (vy(i,j),j=1,nj)
        end do

        write(7,*)  ' radial velocity'
        do i=1,ni
          write(7,103) (0.0,j=1,nj)
        end do

        write(7,*)  ' mach number '
        do i=1,ni
          write(7,103) (vmach(i,j),j=1,nj)
        end do

        write(7,*)  ' static pressure'
        do i=1,ni
          write(7,104) (p(i,j),j=1,nj)
        end do

        write(7,*)  ' density '
        do i=1,ni
          write(7,103) (ro(i,j),j=1,nj)
        end do

        write(7,*)  ' variable plotvar '
        do i=1,ni
          write(7,103) (plotvar(i,j),j=1,nj)
        end do

        write(7,*)  ' del rovx '
        do i=1,ni
          write(7,105)(delrovx(i,j),j=1,nj)
        end do

      end if

      close(7)

      end
