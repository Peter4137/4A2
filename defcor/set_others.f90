      subroutine set_others

      use common_block


      implicit none

! Local stuff
      integer  ::  i, j
!      real, dimension(ni,nj) :: v, tstatic
      real, dimension(i_max,j_max) :: v, tstatic
! This routine calculates secondary flow variables from the primary ones
! at every grid point.

! The primary variables are ro, rovx, rovy and roe

! The secondary variables are the velocity components vx(i,j) and vy(i,j),
! the static pressure p(i,j) and the stagnation enthalpy hstag(i,j).
! Note:  "hstag"  not  "ho".

      do i=1,ni
            do j=1,nj
                  vx(i,j) = rovx(i,j)/ro(i,j)
                  vy(i,j) = rovy(i,j)/ro(i,j)
                  p(i,j) = (gamma-1)*(roe(i,j)-0.5*ro(i,j)*norm2([vx(i,j),vy(i,j)])**2)
                  hstag(i,j) = (roe(i,j)+p(i,j))/ro(i,j)
            end do
      end do
      
! INSERT your code here
      end
