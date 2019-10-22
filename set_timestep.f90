      subroutine set_timestep

      use common_block


      implicit none

! Local stuff
      real  ::  umax

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
      deltat = cfl*dmin/(umax+astag)
! INSERT your code here
      