      subroutine check_grid


      use common_block


      implicit none

! Local stuff
      integer :: i, j,test
      real :: xSum, ySum,small

! Check your grid and areas for both the "bump" and the "bend"
! test data.

! First check that all areas are positive (by program or writing out)

! INSERT your code here

! Next check that the sum of the length vectors of the 4 faces
! of every element is very nearly zero in each coordinate direction.
! It is absolutely essential that this is correct !
! If not go back and check your subroutine "generate_grid".

! Careful with a test of the form
!          if( a == 0.0 ) then .....
! This will probably never be true.  Computers work to a finite number of
! Significant figures and "a" will probably be +0.0000001 or -0.0000001.
! Test for something like
!          if( abs(a) <= small_number ) then ...

! Insert your code here

! Any other tests that you can think of. For example you could plot
! contours of  "area(i,j)" by using  --  call output_hg(area) .

      end
