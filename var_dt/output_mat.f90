      subroutine output_mat(input)

! You do not need to change this subroutine but are welcome to do so.

      use common_block


      implicit none

! Local stuff
      integer ::  i, j, input

! Make up some dummy data for cases when the grid is the only thing of interest.

      if( input==0 ) then
        do i=1,ni
          do j=1,nj
            ro(i,j) = 1.226
            p(i,j)  = 101325.
            vx(i,j) = 0.
            vy(i,j) = 0.
          end do
        end do
      end if

      open(unit=7,file='euler.mat')
      write(7,700) ni,nj
  700 format(i5,1x,i5)
      do i=1,ni
        do j=1,nj
          write(7,701) x(i,j),y(i,j),ro(i,j),vx(i,j),vy(i,j),p(i,j)
        end do
      end do
  701 format(2(1x,f8.5),f8.5,1x,f8.2,1x,f8.2,1x,f10.1)

      close(7)

      end
