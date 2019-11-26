      subroutine set_timestep

      use common_block


      implicit none

! Local stuff
      real  ::  umax
      real :: c_local , v_local, tstag_local
! Tthis subroutine sets the length of the time step based on the
! stagnation speed of sound "astag" and the minimum length scale
! of any element, "dmin". The timestep must be called "deltat"

! An assumption that the maximum flow speed will be equal to "astag"
! is also made. This will be pessimistic for subsonic flows
! but may be optimistic for supersonic flows. In the latter case the
! length of the time step as determined by "cfl" may need to be reduced.

! The cfl number was input as data in data set "flow"

      astag  = sqrt(gamma*rgas*tstagin)
      umax   = astag
!      deltat = cfl*dmin/abs(umax+astag)
! INSERT your code here
      do i = 1,ni-1
            do j = 1, nj-1
                  t_stag_local = (hstag(i,j)+hstag(i+1,j)+hstag(i,j+1)+hstag(i+1,j+1))/4/cp
                  c_local = sqrt(gamma*rgas*(tstag_local))
                  v_local = (norm2([vx(i,j),vy(i,j)]) + &
                              norm2([vx(i+1,j),vy(i+1,j)]) + &
                              norm2([vx(i,j+1),vy(i,j+1)]) + &
                              norm2([vx(i+1,j+1),vy(i+1,j+1)]))/4
                  step(i,j) = cfl*dmin(i,j)/(c_local+v_local)
            enddo
         enddo
      end