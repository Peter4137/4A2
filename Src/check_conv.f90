      subroutine check_conv

! You should not need to change this subroutine

      use common_block


      implicit none

! Local stuff
      integer ::  i, j, imax, jmax
      real    :: delromax, delrovxmax, delrovymax, delroemax, delroavg
      real    :: delrovxavg, delrovyavg, delroeavg, delta, flow_ratio

! This subroutine checks the changes in all primary variables over
! the last 5 steps.
! N.B. the first time this routine is called ro_old, etc. do not
! hold sensible values.  Therefore nothing is written out when nstep = 5.
!
! It also writes a short output summary to the screen and a file for
! plotting the convergence history to units 3 and 31.

      delromax   = 0.0
      delrovxmax = 0.0
      delrovymax = 0.0
      delroemax  = 0.0
      delroavg   = 0.0
      delrovxavg = 0.0
      delrovyavg = 0.0
      delroeavg  = 0.0
      imax = 0
      jmax = 0

! "imax,jmax" is the grid point where the change in rovx is a max.

      do i=1,ni
        do j=1,nj

          delta = abs(ro(i,j) - ro_old(i,j))
          if(delta > delromax) delromax = delta
          delroavg = delroavg + delta

          delta = abs(rovx(i,j)-rovx_old(i,j))
          if(delta > delrovxmax) then
            delrovxmax = delta
            imax = i
            jmax = j
          end if
          diffrovx(i,j) = delta
          delrovxavg = delrovxavg + delta

          delta = abs(rovy(i,j) - rovy_old(i,j))
          if(delta > delrovymax) delrovymax = delta
          delrovyavg = delrovyavg + delta

          delta = abs(roe(i,j) - roe_old(i,j))
          if(delta > delroemax) delroemax = delta
          delroeavg = delroeavg + delta

        end do
      end do

! Calculate the average changes

      delroavg   =  delroavg/ncells/ref_ro
      delrovxavg = delrovxavg/ncells/ref_rovx
      delrovyavg = delrovyavg/ncells/ref_rovy
      delroeavg  = delroeavg/ncells/ref_roe
      delrovxmax = delrovxmax/ref_rovx
      delrovymax = delrovymax/ref_rovy

      emax = amax1(delrovxmax,delrovymax)
      eavg = amax1(delrovxavg,delrovyavg)

! Store the maximum change in rovx as emax to be printed out.

! Save the current values of the primary variables as prop_old values
! for use in the next convergence check.

      do i=1,ni
        do j=1,nj
          ro_old(i,j)   = ro(i,j)
          rovx_old(i,j) = rovx(i,j)
          rovy_old(i,j) = rovy(i,j)
          roe_old(i,j)  = roe(i,j)
        end do
      end do

! Write the average changes in the primary variables to units 3 and 31
! for use in convergence plotting (pltconv and paraview respectively).

      if( nstep > 5) then
        write(3,300) delroavg, delrovxavg, delrovyavg, delroeavg
  300   format(4e13.6)
        write(31,"(i5,a1,1x,4(f13.6,a1,1x))")  &
             nstep,',',delroavg,',',delrovxavg,',',delrovyavg,',',  &
             delroeavg
      end if

! Write a short output summary to the screen.

      flow_ratio = flow(ni)/flow(1)
      write(*,*) ' time step number ', nstep
      write(*,600) emax,imax,jmax,eavg
  600 format(' emax= ',e10.3,' at imax = ',i5,' jmax= ',i5,' eavg= ', e10.3)
      write(6,*) 'inlet flow= ',flow(1),' outlet to inlet flow ratio', flow_ratio
      write(6,*) 'outlet flow= ',flow(ni)

      end
