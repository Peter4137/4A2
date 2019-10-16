      subroutine generate_grid

      use common_block


      implicit none

! Local variables
      integer :: i, j

! Calculate x and y values  x(i,j),y(i,j) of the grid nodes.

! For each value of "i" the i-grid line joins (xlow(i),ylow(i)) to
! (xhigh(i),yhigh(i)). for each value of "i" grid points (nodes) should be
! linearly interpolated between these values for j between 1 and nj.
! i.e.  x(i,1) should be xlow(i), x(i,nj) should be xhigh(i), etc.

      REAL, DIMENSION(1:ni,nj) :: x, y

      do i=1,ni 
            do j=0,nj-1
                  x(i,j) = xlow(i)
                  y(i,j) = ylow(i) + j*(yhigh(i)-ylow(i))/nj
            end do
      end do

! INSERT your code here

! Calculate the areas of the cells area(i,j)
! (N.B. there are (ni-1) x (nj-1) cells.

! The area of a quadrilateral (regular or irregular) can be shown to be
! half of the cross product of the vectors forming the diagonals.
! see Hirsch volume 1, section 6.2.1. (or lecture).
! Make sure that the area comes out positive!
      REAL, DIMENSION(1:ni-1,nj-1) :: area
      do i=1,ni-1
            do j=1,nj-1
                  area(i,j) = 0.5*(((x(i+1,j+1)-x(i,j))*(y(i,j+1)-y(i+1,j)))-((y(i+1,j+1)-y(i,j))*(x(i,j+1)-x(i+1,j))))
                  if(area(i,j).lt.0) then
                        area(i,j) = -area(i,j)
                  end if
            end do
      end do
! INSERT your code here

! Calculate the x and y components of the length vector of the i-faces
! (i.e. those corresponding to i = constant).
! The length vector of a face is a vector normal to the face wi
! magnitude equal to the length of the face.
! It is positive in the direction of an inward normal to the cell i,j .
! Call these lengths dlix(i,j) and dliy(i,j)
      REAL, DIMENSION(1:ni-1,nj-1) :: dlix,dliy,dljx,dljy
      REAL :: dmin
      dmin = 1000
      do i=1,ni-1
            do j=1,nj-1
                  dlix(i,j) = y(i,j+1)-y(i,j)

                  if(dlix(i,j).lt.dmin) then
                        dmin = dlix(i,j)
                  end if

                  dliy(i,j) = x(i,j+1)-x(i,j)

                  if(dliy(i,j).lt.dmin) then
                        dmin = dliy(i,j)
                  end if

                  dljx(i,j) = y(i+1,j)-y(i,j)

                  if(dljx(i,j).lt.dmin) then
                        dmin = dljx(i,j)
                  end if

                  dljy(i,j) = x(i+1,j)-x(i,j)

                  if(dljy(i,j).lt.dmin) then
                        dmin = dljy(i,j)
                  end if

            end do
      end do

! INSERT your code here

! Now calculate the x and y components of the length vector of the j-faces. (i.e. those corresponding to j = constant)
! Call these lengths dljx(i,j) and dljy(i,j)

! INSERT your code here

! Now find "dmin" the minimum length scale of any element. This is
! defined as the length of the shortest side of the element.
! Call this minimum "dmin". it is used to set the time step from the cfl no.

! Insert your code here (or in the do loops above).

      write(6,*)  ' overall minimum element size = ', dmin

      end
