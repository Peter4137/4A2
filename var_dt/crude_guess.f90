      subroutine crude_guess

! You should not need to touch this subroutine

      use common_block

      implicit none

! This subroutine makes a very simple and inaccurate guess of the
! flow field. enough to get the program working.

! Local stuff
      integer i, j, jmid
      real    tdown, vdown, rodown, dx, dy, ds, xvel, yvel

      jmid   = nj/2
      tdown  = tstagin*(pdown/pstagin)**fga
      vdown  = sqrt(2*cp*(tstagin - tdown))
      rodown = pdown/rgas/tdown

      do j=1, nj
        do i=1, ni-1
           dx  = x(i+1,jmid) - x(i,jmid)
           dy  = y(i+1,jmid) - y(i,jmid)
           ds  = sqrt(dx*dx + dy*dy)
           xvel      = vdown*dx/ds
           yvel      = vdown*dy/ds
           rovx(i,j) = rodown*xvel
           rovy(i,j) = rodown*yvel
           ro(i,j)   = rodown
           roe(i,j)  = rodown*(cv*tdown + 0.5*vdown*vdown)
        end do

        rovx(ni,j) = rovx(ni-1,j)
        rovy(ni,j) = rovy(ni-1,j)
        ro(ni,j)   = ro(ni-1,j)
        roe(ni,j)  = roe(ni-1,j)

      end do

      end
