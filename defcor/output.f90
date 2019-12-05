      subroutine output(input)

! You should not need to change this subroutine.  Paraview works with 3d
! flowfields.

      use common_block


      implicit none

! Local stuff
      integer  :: i, j, k, input
      real    :: q(5,i_max,j_max,2), xs(3,i_max,j_max,2), kenergy

! Set up the coordinates, and insert a fictitious 3rd dimension

      do i=1,ni
        do j=1,nj
          xs(1,i,j,1) = x(i,j)
          xs(2,i,j,1) = y(i,j)
          xs(3,i,j,1) = 0.0

          xs(1,i,j,2) = x(i,j)
          xs(2,i,j,2) = y(i,j)
          xs(3,i,j,2) = 0.1
        end do
      end do

! Take the input switch (1 for flow data, 0 for fake data)
! Fake data is necessary in the early stages of development to check grids, when
! there may not be a solution

      if(input == 1) then

        do k=1,2
          do j=1,nj
            do i=1,ni

              q(1,i,j,k) = ro(i,j)
              q(2,i,j,k) = rovx(i,j) / ro(i,j)
              q(3,i,j,k) = rovy(i,j) / ro(i,j)
              q(4,i,j,k) = 0.0

              kenergy    = 0.5 * ro(i,j) * &
                 (q(2,i,j,k)**2 + q(3,i,j,k)**2 + q(4,i,j,k)**2 )

              q(5,i,j,k) = (gamma-1) * (roe(i,j) - kenergy)

            end do
          end do
        end do

      else if(input == 0) then

        do k=1,2
          do j=1,nj
            do i=1,ni

              q(1,i,j,k) = 1.226
              q(2,i,j,k) = 0.0
              q(3,i,j,k) = 0.0
              q(4,i,j,k) = 0.0
              q(5,i,j,k) = 101300.0

            end do
          end do
        end do

      end if

! Open the file for outputting csv data

      open(unit=7,file='euler.csv')

! Write the csv data for paraview

      write(7,*) ' x,y,z,rho,u,v,w,p'
      do k = 1, 2
        do j = 1, nj
          do i = 1, ni
            write(7,"(8(f13.6,a1,1x))") &
              xs(1,i,j,k),',', xs(2,i,j,k),',', xs(3,i,j,k),',', &
              q(1,i,j,k),',', q(2,i,j,k),',', q(3,i,j,k),',',    &
              q(4,i,j,k),',', q(5,i,j,k)
          end do
        end do
      end do

!     close the unit

      close(7)

      end
