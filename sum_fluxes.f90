      subroutine sum_fluxes(iflux,jflux,delprop,prop_inc)

      use common_block


      implicit none

! This subroutine sums the fluxes for each element, calculates the changes
! in the variable "prop" appropriate to a cell (delprop) and distributes them to the
! four corners of the element (stored in prop_inc).

      real, dimension(i_max,j_max) ::  iflux, jflux, prop_inc, delprop,  &
                                       previous, store
!      real :: frkut

! Local stuff
      integer ::    i, j

! Note previous changes prior to updating them (used with second order timesteps)
! N.B. These are cell quantities.  N.B. There are ni-1 x nj-1 cells.

      do i=1,ni-1
        do j=1,nj-1
          previous(i,j) = delprop(i,j)
        end do
      end do

! Find the change in the variable "prop" in each cell over the
! time step "deltat" and save it in "delprop(i,j)".
      do i=1,ni-1
            do j=1,nj-1
                  delprop(i,j) = (step(i,j)/area(i,j))*(iflux(i,j) - iflux(i+1,j) + jflux(i,j) - jflux(i,j+1))
            end do
      end do
! INSERT your code here to calculate the change in the variable
! "prop" over the time step "deltat" and set the result to "delprop(i,j)"


! Now distribute the changes equally to the four corners of each
! cell (nodes). Each interior grid point receives one quarter of the change
! from each of the four cells adjacent to it.  Edge nodes do not have four adjacent
! cells.

      do i=2,ni-1
        do j=2,nj-1
! INSERT your code here to calculate "prop_inc" the change to be added to
! each interior node.
          prop_inc(i,j) = 0.25*(delprop(i,j)+delprop(i-1,j)+delprop(i,j-1)+delprop(i-1,j-1))
        enddo
      enddo

! Now deal with the changes to the upper and lower boundaries.
! These receive half the change from each of the two cells adjacent to them.

      do i=2,ni-1

! INSERT your code here to calculate "prop_inc" for the nodes with j=1.
        prop_inc(i,1) = 0.5*(delprop(i,1)+delprop(i-1,1))

! INSERT your code here to calculate "prop_inc" for the nodes with j=nj.
        prop_inc(i,nj) = 0.5*(delprop(i,nj-1)+delprop(i-1,nj-1))

      end do

! Now deal with changes to the inlet & outlet boundary points.
! These receive half the change from each of the two cells adjacent to them.

      do j=2,nj-1

! INSERT your code here to calculate "prop_inc" for the nodes with i=ni.
        prop_inc(ni,j) = 0.5*(delprop(ni-1,j)+delprop(ni-1,j-1))

! INSERT your code here to calculate "prop_inc" for the nodes with i=1.
        prop_inc(1,j) = 0.5*(delprop(1,j)+delprop(1,j-1))

      end do

! Finally find the changes to be added to the four corner points.
! These receive the full change from the single cell of which they form one corner.

! INSERT your code here to calculate "prop_inc" for the node with i=1,j=1.
     prop_inc(1,1) = delprop(1,1)

! INSERT your code here to calculate "prop_inc" for the node with i=1,j=nj.
     prop_inc(1,nj) = delprop(1,nj-1)

! INSERT your code here to calculate "prop_inc" for the node with i=ni,j=1.
     prop_inc(ni,1) = delprop(ni-1,1)

! INSERT your code here to calculate "prop_inc" for the node with i=ni,j=nj.
     prop_inc(ni,nj) = delprop(ni-1,nj-1)

! If the values in delprop have been hi-jacked by the second order timestep extension,
! restore them here to the true residuals. (The second order extension should have
! loaded them into the variable "store" bwfore hi-jacking them)
! These are used in the convergence check and also in future improvements to the scheme.

      end
