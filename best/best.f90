      subroutine apply_bconds

      use common_block

      implicit none

! Local stuff
      integer ::  j
      real ::  rfin, rfin1, rostagin, gm1
      real ::  tstat, eke, vel

! This subroutine applies the boundary conditions that p = pdown
! at i = ni. At the inlet boundary the change in density is relaxed
! to obtain roinlet(j) which is then used to obtain the other properties
! at inlet assuming isentropic flow from stagnation conditions "pstagin"
! and tstagin" together with the specified inlet flow angle "alpha1".

! Because the inlet condition may become unstable, it is safer to
! relax the changes in inlet density by a factor "rfin"
! Typically "rfin" = 0.25 as set below. Reduce this if the inlet
! becomes unstable.

! It is also worth checking if "roinlet" is greater than "rostagin"
! and setting roinlet to 0.9999*rostagin if it is.
! This saves the program crashing during severe transients.

      rfin     = 0.25
      rfin1    = 1.0-rfin
      rostagin = pstagin/(rgas*tstagin)
      gm1      = gamma - 1.0

      do j=1,nj

        if (nstep==1) then
          roinlet(j) = ro(1,j)
        else
          roinlet(j) = rfin*ro(1,j) + rfin1*roinlet(j)
        endif

        if ( roinlet(j) > 0.9999*rostagin ) then
          roinlet(j) = 0.9999*rostagin
        endif
        
        tstat = tstagin*(roinlet(j)/rostagin)**gm1

        vel = sqrt(2*cp*(tstagin-tstat))
        
        rovx(1,j) = roinlet(j)*vel*cos(alpha1)
        rovy(1,j) = roinlet(j)*vel*sin(alpha1)
        
        roe(1,j) = roinlet(j)*(cv*tstat + 0.5*vel**2)
        p(1,j) = gm1*(roe(1,j) - 0.5*roinlet(j)*vel**2)
        vx(1,j) = rovx(1,j)/roinlet(j)
        vy(1,j) = rovy(1,j)/roinlet(j)
        hstag(1,j) = cp*tstat + 0.5*vel**2

! INSERT your code here to calculate p(1,j), rovx(1,j), rovy(1,j)
! and roe(1,j)  from roinlet(j), pstagin, tstagin  and alpha1.
! Also set vx(1,j), vy(1,j) and hstag(1,j)

      end do

! Set the pressure at the downstream boundary i = ni to the exit
! static pressure "pdown" for all j values.

      do j=1,nj
            p(ni,j) = pdown
      end do

      ! INSERT your code here

      end
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

      end
      subroutine check_grid


      use common_block


      implicit none

! Local stuff
      integer :: i, j,test
      real :: xSum, ySum,small

! Check your grid and areas for both the "bump" and the "bend"
! test data.

! First check that all areas are positive (by program or writing out)

      do i=1,ni-1
            do j=1,nj-1
                  if(area(i,j).lt.0) then
                        write(6,*) 'Areas are not all positive'
                        exit
                  endif
            end do
      end do
                  
! INSERT your code here

! Next check that the sum of the length vectors of the 4 faces
! of every element is very nearly zero in each coordinate direction.
! It is absolutely essential that this is correct !
! If not go back and check your subroutine "generate_grid".

      small = 0.000002
      do i=1,ni-1
            do j=1,nj-1
                  xSum = dlix(i,j) + dljx(i,j) - dlix(i+1,j) - dljx(i,j+1)
                  if (abs(xSum).GT.small) then
                        write(6,*) 'X-grid not zero!'
                        write(6,*) xSum
                        stop
                  endif
                  xSum = dliy(i,j) + dljy(i,j) - dliy(i+1,j) - dljy(i,j+1)
                  if (abs(ySum).GT.small) then
                        write(6,*) 'Y-grid not zero!'
                        write(6,*) ySum
                        stop
                  endif
            end do
      end do

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

module common_block

! Declare all variables and arrays for use in main and subroutines.
! The variables are grouped according to type of use.  Some are needed only
! for specific extensions.

      integer, parameter    :: i_max=201, j_max=51

! Variables associated with the grid
      real, dimension(i_max) :: xlow, ylow, xhigh, yhigh
      real, dimension(i_max,j_max) :: x, y, area, dli, dlj, dljx, dljy, dlix, dliy
      real, dimension(i_max,j_max)  ::  dmin
      integer ::   ni,nj
      real, dimension(i_max,j_max) :: step

! Variables to hold node increments
      real, dimension(i_max,j_max) ::  ro_inc, roe_inc, rovx_inc, rovy_inc

! Variables to implement Runge Kutta option
      real, dimension(i_max,j_max) ::  ro_start, roe_start, rovx_start, rovy_start
      real :: frkut
      integer  ::  nrkut, nrkut_max

! Title of test case (read in with data)
      character(LEN=80) :: title

! Gas Properties
      real ::    rgas, gamma, cp, cv, fga

! Primary Variables (Values at nodes)
      real, dimension(i_max,j_max)  ::  ro, rovx, rovy, roe

! Residuals (changes in cells)
      real, dimension(i_max,j_max)  ::  delro, delrovx, delrovy, delroe

! Convergence check
      real, dimension(i_max,j_max)  ::  ro_old, rovx_old, rovy_old, roe_old, diffrovx

! Smoothing variables for deferred correction option
      real, dimension(i_max,j_max)  ::  corr_ro, corr_rovx, corr_rovy, corr_roe
      real  ::  fcorr

! Other variables
      real, dimension(i_max,j_max)  ::  p, hstag, vx, vy

! Fluxes across cell faces
      real, dimension(i_max,j_max) ::   fluxi_mass, fluxj_mass, fluxi_xmom, fluxj_xmom, &
                                        fluxi_ymom, fluxj_ymom, fluxi_enth, fluxj_enth
      real, dimension(i_max) :: flow

! Timestep and other stuff
      real ::  cfl, smooth_fac, deltat, conlim, emax,  &
               eavg, grid_ratio, facsec
      integer ::  nstep, nsteps, ncells

! Boundary conditions
      real, dimension(j_max) :: pin, roinlet
      real  ::  pstagin, tstagin, alpha1, astag, pdown

! Useful reference values
      real ::  ref_p, ref_ro, ref_t, ref_v, ref_rovx, ref_rovy, ref_roe, roin

! Initial guess variables
      real, dimension(i_max) ::  aflow, v_guess, ro_guess

! Local stuff
!      integer i,j

end module common_block
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
      program euler

      use common_block

      implicit none



! Local stuff
      integer i,j

! Open files to store the convergence history. Plotting is done via a separate
! program. "euler.log" is for use by pltconv. "pltcv.csv" is for use by paraview.

      open(unit=3,file='euler.log')
      open(unit=31,file='pltcv.csv')

! "read_data": to read in the data on the duct and geometry and flow conditions.

      call read_data

! "generate_grid": to set up the grid coordinates, element areas and projected
! lengths of the sides of the elements.

      call generate_grid

! You can call subroutines "output*" here to plot out the grid you have generated.
! "output" writes "euler.csv" for paraview, "ouput_hg" write "euler.plt" for eulplt,
! "output_mat" writes "euler.mat" for matlab.

      call output(0)
      call output_hg(p,0)
      call output_mat(0)

! "check_grid": to check that the areas and projected lengths are correct.

      call check_grid

! "crude_guess" is what its name says. it enables you to
! start a calculation and obtain a solution but it will take longer than
! necessary. when your program is working you should replace it
! with "flow_guess" to obtain a better guess and a faster solution.

      ! call crude_guess
      ! call flow_guess
      call new_guess
      ! stop

! You can call "output" here to plot out your initial guess of
! the flow field.

      call output(1)
      call output_mat(1)

! "set_timestep": to set the length of the timestep.
! initially this is a constant time step based on a conservative guess
! of the mach number.

      call set_timestep

!************************************************************************
!     start the time stepping do loop for "nsteps" loops.
!************************************************************************
    
      do nstep = 1, nsteps

            do i=1,ni
              do j=1,nj
                ro_start(i,j) = ro(i,j)
                roe_start(i,j) = roe(i,j)
                rovx_start(i,j) = rovx(i,j)
                rovy_start(i,j) = rovy(i,j)
              end do
            end do
       
    ! Runge kutta scheme
            nrkut_max = 4
            do nrkut = 1,nrkut_max
                frkut = 1.0/(1+nrkut_max-nrkut)
    
    ! "set_others" to set secondary flow variables.
                
                call set_others
    
    ! "apply_bconds" to apply inlet and outlet values at the boundaries of the domain.
                
                call apply_bconds
    
    ! "set_fluxes" to set the fluxes of the mass, momentum, and energy throughout the domain.
                
                call set_fluxes
    
    ! "sum_fluxes" applies a control volume analysis to enforce the finite volume method
    ! for each cell (calculates the residuals) and sets the increments for the nodal values.
                
                call sum_fluxes(fluxi_mass,fluxj_mass,delro  , ro_inc, frkut)
                call sum_fluxes(fluxi_enth,fluxj_enth,delroe ,roe_inc, frkut)
                call sum_fluxes(fluxi_xmom,fluxj_xmom,delrovx,rovx_inc, frkut)
                call sum_fluxes(fluxi_ymom,fluxj_ymom,delrovy,rovy_inc, frkut)
    !
    ! Update solution
    
                do i=1,ni
                do j=1,nj
                      ro  (i,j) = ro_start  (i,j) + ro_inc  (i,j)
                      roe (i,j) = roe_start (i,j) + roe_inc (i,j)
                      rovx(i,j) = rovx_start(i,j) + rovx_inc(i,j)
                      rovy(i,j) = rovy_start(i,j) + rovy_inc(i,j)
                end do
                end do
    
    ! Smooth the problem to ensure it remains stable.
                
                call smooth(ro, corr_ro)
                call smooth(rovx, corr_rovx)
                call smooth(rovy, corr_rovy)
                call smooth(roe, corr_roe)
            end do 
    ! Check convergence and write out summary every 5 steps
     
        if(mod(nstep,5)==0) then
            call set_timestep
            call check_conv
        end if

! Stop looping if converged to the input tolerance "conlim"

        if( emax < conlim .and.  eavg < (0.5*conlim) ) then
          write(6,*) ' Calculation converged in ',nstep,' iterations'
          write(6,*) ' To a convergence limit of ', conlim
          exit
        endif
      end do
      !call output(1)
      !call output_mat(1)

!************************************************************************
!  end of time stepping do loop for "nsteps" loops.
!************************************************************************

! Calculation finished. call "output" to write the plotting file.
! N.B. Solution hasn't necessarily converged.

      call output(1)
      call output_hg(p,1)
      call output_mat(1)
!
      close(3)
      close(31)


      end
      subroutine flow_guess


      use common_block


      implicit none

! Local stuff
      real ::   mflow, machlim, tlim, mach_num, Tdown,e,dy,dx,dxy, rodown, vdown, gm1
      real, dimension(i_max) :: Tstatic,Pstatic
      integer ::  i, j

      gm1 = gamma-1
! In this subroutine we make an initial guess of the primary variables
! i.e.  ro, rovx, rovy and roe.
! The guess does not need to be very accurate but the better
! it is the faster your program will converge.
! You should assign values to ro(i,j), rovx(i,j), rovy(i,j) and roe(i,j)
! at every grid point in this subroutine.

! Work out the length of each "i" line between grid points "i,1" and "i,nj"
! and call it  "aflow(i)" .
      do i=1,ni
            aflow(i) = norm2([x(i,nj)-x(i,1), y(i,nj)-y(i,1)])
      end do

! INSERT your code here

! Make an initial guess of the density and velocity at the exit by
! assuming isentropic flow between the inlet stagnation pressure pstagin
! and temperature tstagin and the exit static pressure pdown.
! Use these together with "aflow(ni)" to estimate the mass flow rate.
! call this "mflow".
      Tdown = tstagin*(pdown/pstagin)**(gm1/gamma)
      rodown = (pdown/(rgas*Tdown))
      vdown = sqrt(2*cp*(tstagin-Tdown))
      mflow = rodown*aflow(ni)*vdown
! INSERT your code here

! Set a limit to the maximum allowable mach number in the initial
! guess. call this "machlim". calculate the corresponding temperature.

      machlim = 1.0
      tlim = tstagin/(1.0 + 0.5*(gamma-1.0)*machlim*machlim)

! Now estimate the velocity and density at every "i" line.
! Call the velocity v_guess(i) and the density ro_guess(i).


! First assume that the density is constant and equal to the exit
! density calculated above and that the flow is perpendicular to the
! "i" = constant lines and hence occupies artea aflow(i).
! Use continuity to estimate the flow velocity v_guess(i).
! Use this velocity to calculate the static temperature assuming
! that the stagnation temperature is constant.
! Check that this temperature is not less than tlim and set = tlim
! if it is.
! Next use this temperature and isentropic flow to obtain a better
! estimate of the density, ro_guess(i).
! Use this density and continuity to obtain a better estimate of
! the velocity, set = v_guess(i).

      do i=1,ni
            v_guess(i) = mflow/(rodown*aflow(i))
            Tstatic(i) = tstagin - (v_guess(i)**2)/(2*cp)
            if(Tstatic(i).lt.tlim) then
                  Tstatic(i) = tlim
            end if
            ro_guess(i) = (pstagin*(Tstatic(i)/tstagin)**(gamma/gm1))/(rgas*Tstatic(i))
            v_guess(i) = mflow/(ro_guess(i)*aflow(i))
      end do

! INSERT your code here

! Direct the velocity found above along the j= constant grid lines to find
! the velocity vx(i,j) in the  x  direction and vy(i,j) in the y.
! Use these and ro_guess(i) to set rovx(i,j), rovy(i,j) and roe(i,j).
! Also set ro(i,j).
! Note that roe(i,j) includes the kinetic energy component of the
! internal energy.

      do i=1,ni                                             
            do j=1,nj
                  dx  = x(i,nj)-x(i,1)
                  dy  = y(i,nj)-y(i,1)
                  dxy = norm2([dx,dy])
                  vx(i,j) = v_guess(i)*dy/dxy
                  vy(i,j) = -v_guess(i)*dx/dxy
                  ro(i,j) = ro_guess(i)
                  rovx(i,j) = ro_guess(i)*vx(i,j) 
                  rovy(i,j) = ro_guess(i)*vy(i,j)
                  roe(i,j) = ro_guess(i)*(cv*Tstatic(i)+0.5*v_guess(i)**2)
            end do
      end do
      


! INSERT your code here
     
! Store the "old" values of the variables for use in the first
! convergence check in subroutine "check_conv"

      do i=1,ni
        do j=1,nj
          ro_old(i,j)   = ro(i,j)
          rovx_old(i,j) = rovx(i,j)
          rovy_old(i,j) = rovy(i,j)
          roe_old(i,j)  = roe(i,j)
        end do
      end do

      end
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

      ! REAL, DIMENSION(1:ni,nj) :: x, y

      do i=1,ni 
            do j=1,nj
                  ! x(i,j) = ((xhigh(i)-xlow(i))/(nj-1))*(j-1) + xlow(i)
                  ! y(i,j) = ((yhigh(i)-ylow(i))/(nj-1))*(j-1) + ylow(i)
                  x(i,j) = xlow(i) + j*(xhigh(i)-xlow(i))/nj
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
      ! REAL, DIMENSION(1:ni-1,nj-1) :: area
      do i=1,ni-1
            do j=1,nj-1
                  area(i,j) = abs(0.5*(((x(i+1,j+1)-x(i,j))*(y(i,j+1)-y(i+1,j)))-((y(i+1,j+1)-y(i,j))*(x(i,j+1)-x(i+1,j)))))
            end do
      end do
! INSERT your code here

! Calculate the x and y components of the length vector of the i-faces
! (i.e. those corresponding to i = constant).
! The length vector of a face is a vector normal to the face wi
! magnitude equal to the length of the face.
! It is positive in the direction of an inward normal to the cell i,j .
! Call these lengths dlix(i,j) and dliy(i,j)
      ! REAL, DIMENSION(1:ni-1,nj-1) :: dlix,dliy,dljx,dljy
      ! REAL :: dmin
      do i=1,ni
            do j=1,nj-1
                  dlix(i,j) = y(i,j+1)-y(i,j)
                  dliy(i,j) = x(i,j)-x(i,j+1)

                  dmin(i,j) = norm2([dlix(i,j), dliy(i,j)])
                  
            end do
      end do

      do i=1,ni-1
            do j=1,nj
                  dljx(i,j) = y(i,j)-y(i+1,j)
                  dljy(i,j) = x(i+1,j)-x(i,j)
                  if(dmin(i,j)>0.0001) then
                        dmin(i,j) = min(dmin(i,j), norm2([dljx(i,j), dljy(i,j)]))
                  else
                        dmin(i,j) = norm2([dljx(i,j), dljy(i,j)])
                  endif 
                  if(dmin(i,j).lt.0.00001)then
                        write(6,*) dlix(i,j), dliy(i,j), i, j
                  endif
            end do
      end do

      ! write(6,*) 

! INSERT your code here

! Now calculate the x and y components of the length vector of the j-faces. (i.e. those corresponding to j = constant)
! Call these lengths dljx(i,j) and dljy(i,j)

! INSERT your code here

! Now find "dmin" the minimum length scale of any element. This is
! defined as the length of the shortest side of the element.
! Call this minimum "dmin". it is used to set the time step from the cfl no.

! Insert your code here (or in the do loops above).

      !write(6,*)  ' overall minimum element size = ', dmin

      end
    subroutine new_guess

    use common_block
    !
    !     This subroutine makes a guess of the flow field assuming no
    !     variation in the cross-flow (j) direction and assuming a
    !     linear variation of flow properties from inlet to outlet.
    !
    implicit none

    real :: ds, dx, dy, pinlet, rodown, rolocal, tdown
    real :: tinlet, tlocal, vdown, vinlet, vlocal, xvel, yvel

    integer :: i,j,jmid

    jmid   = nj/2
    tdown  = tstagin*(pdown/pstagin)**fga
    vdown  = sqrt(2*cp*(tstagin - tdown))
    rodown = pdown/rgas/tdown
    pinlet = 55000.
    tinlet  = tstagin*(pinlet/pstagin)**fga
    vinlet  = sqrt(2*cp*(tstagin - tinlet))
    roin    = pinlet/rgas/tinlet
    write(6,*) sqrt(gamma*tinlet*rgas)
    do j=1,nj
    do i=1,ni-1
        dx  = x(i+1,jmid) - x(i,jmid)
        dy  = y(i+1,jmid) - y(i,jmid)
        ds  = sqrt(dx*dx + dy*dy)
        vlocal    = vinlet  + (vdown-vinlet)*float(i-1)/float(ni-1)
        rolocal   = roin    + (rodown-roin)*float(i-1)/float(ni-1)
        tlocal    = tinlet  + (tdown-tinlet)*float(i-1)/float(ni-1)
        
        vx(i,j)      = vlocal*dx/ds
        vy(i,j)      = vlocal*dy/ds
        rovx(i,j) = rolocal*vx(i,j)
        rovy(i,j) = rolocal*vy(i,j)
        ro(i,j)   = rolocal
        roe(i,j)  = rolocal*(cv*tlocal + 0.5*vlocal*vlocal)
        !if(j==25)then
        !        write(6,*) norm2([vx(i,j), vy(i,j)]) 
        !end if
        end do
    rovx(ni,j) = rovx(ni-1,j)
    rovy(ni,j) = rovy(ni-1,j)
    ro(ni,j)  = ro(ni-1,j)
    roe(ni,j) = roe(ni-1,j)
    end do
    
    do i=1,ni
        do j=1,nj
          ro_old(i,j)   = ro(i,j)
          rovx_old(i,j) = rovx(i,j)
          rovy_old(i,j) = rovy(i,j)
          roe_old(i,j)  = roe(i,j)
        end do
      end do

    end subroutine new_guess
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
      subroutine output_hg(plotvar,input)

! You should not need to touch this subroutine.

! This subroutine writes a file to unit 7 for use by the plotting
! program "eulplt".

      use common_block


      implicit none

! Local stuff
      integer ::   input, i,j

      real, dimension(i_max,j_max)  ::   pstag, vmach, plotvar
      real       eke,tstat,velsq,tstag

      open(unit=7,file='euler.plt')

      write(7,100) title
100   format(a80)

      cp = rgas*gamma/(gamma-1.)
      cv = cp/gamma
      write(7,101) 1,ni,nj,0,0,1,ni,0,0
      write(7,102) cp,gamma

      do i=1,ni
        write(7,102) (x(i,j),j=1,nj)
        write(7,102) (y(i,j),j=1,nj)
      end do
101   format(16i5)
102   format(10f10.5)

! Calculate the secondary variables

      if( input==1 ) then
        do i=1,ni
          do j=1,nj
            vx(i,j)  = rovx(i,j)/ro(i,j)
            vy(i,j)  = rovy(i,j)/ro(i,j)
            eke      = 0.5*(vx(i,j)*vx(i,j) + vy(i,j)*vy(i,j))
            tstat    = (hstag(i,j) - eke)/cp
            p(i,j)   = ro(i,j)*rgas*tstat
            roe(i,j) = ro(i,j)*(cv*tstat + eke)
          end do
        end do

! Calculate the mach number and stagnation presssure

        do i=1,ni
          do j=1,nj
            tstat = p(i,j)/rgas/ro(i,j)
            velsq = vx(i,j)*vx(i,j) + vy(i,j)*vy(i,j)
            tstag = tstat + 0.5*velsq/cp
            vmach(i,j) = sqrt(velsq/(gamma*rgas*tstat))
            pstag(i,j) = p(i,j)*(tstag/tstat)**(1/fga)
          end do
        end do

        write(7,*) ' time step number', nstep

        write(7,*)  ' axial velocity'
        do i=1,ni
          write(7,103) (vx(i,j),j=1,nj)
        end do
103     format(10f10.4)

        write(7,*)  ' y  velocity'
        do i=1,ni
          write(7,103) (vy(i,j),j=1,nj)
        end do

        write(7,*)  ' radial velocity'
        do i=1,ni
          write(7,103) (0.0,j=1,nj)
        end do

        write(7,*)  ' mach number '
        do i=1,ni
          write(7,103) (vmach(i,j),j=1,nj)
        end do

        write(7,*)  ' static pressure'
        do i=1,ni
          write(7,104) (p(i,j),j=1,nj)
        end do
104     format(10f10.1)

        write(7,*)  ' density '
        do i=1,ni
          write(7,103) (ro(i,j),j=1,nj)
        end do

        write(7,*)  ' variable plotvar '
        do i=1,ni
          write(7,103) (plotvar(i,j),j=1,nj)
        end do

        write(7,*)  ' del rovx '
        do i=1,ni
          write(7,105) (delrovx(i,j),j=1,nj)
        end do
105     format(10e10.4)

! This is the dummy part used only for initial stages of grid debugging
      else

        nstep = -1
        do i=1,ni
          do j=1,nj
            ro(i,j)  = 1.226
            vx(i,j)  = 0.
            vy(i,j)  = 0.
            tstat    = 288.15
            p(i,j)   = ro(i,j)*rgas*tstat
            roe(i,j) = ro(i,j)*cv*tstat
          end do
        end do

! Calculate the mach number and stagnation presssure

        do i=1,ni
          do j=1,nj
            tstat = p(i,j)/rgas/ro(i,j)
            velsq = vx(i,j)*vx(i,j) + vy(i,j)*vy(i,j)
            tstag = tstat + 0.5*velsq/cp
            vmach(i,j) = sqrt(velsq/(gamma*rgas*tstat))
            pstag(i,j) = p(i,j)*(tstag/tstat)**(1/fga)
          end do
        end do

        write(7,*) ' time step number', nstep

        write(7,*)  ' axial velocity'
        do i=1,ni
          write(7,103) (vx(i,j),j=1,nj)
        end do

        write(7,*)  ' y  velocity'
        do i=1,ni
          write(7,103) (vy(i,j),j=1,nj)
        end do

        write(7,*)  ' radial velocity'
        do i=1,ni
          write(7,103) (0.0,j=1,nj)
        end do

        write(7,*)  ' mach number '
        do i=1,ni
          write(7,103) (vmach(i,j),j=1,nj)
        end do

        write(7,*)  ' static pressure'
        do i=1,ni
          write(7,104) (p(i,j),j=1,nj)
        end do

        write(7,*)  ' density '
        do i=1,ni
          write(7,103) (ro(i,j),j=1,nj)
        end do

        write(7,*)  ' variable plotvar '
        do i=1,ni
          write(7,103) (plotvar(i,j),j=1,nj)
        end do

        write(7,*)  ' del rovx '
        do i=1,ni
          write(7,105)(delrovx(i,j),j=1,nj)
        end do

      end if

      close(7)

      end
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
      subroutine read_data

      use common_block


      implicit none

! Local stuff
      real ::   conlim_in, smooth_fac_in
      integer :: i, jmid

! Assign unit 1 to the file 'geom'
! Assign unit 2 to the file 'flow'

      open( unit = 1, file = "../test5_geom")
      open( unit = 2, file = "../test5_flow")
! INSERT your code here

! Read in the title and ni and nj from unit 1.
! Before you do this, check the format of the data file "geom" and read the
! appropriate values for ni & nj.

      read(1,*) title
      read(1,*) NI, NJ

! INSERT your code here

! Check that ni and nj are less than the array dimensions i_max and j_max.
! Then read in xlow, ylow, xhigh, yhigh for each i between 1 and ni.

      ! REAL, DIMENSION(1:ni) :: xlow, ylow, xhigh, yhigh

      if(ni.le.i_max.and.nj.le.j_max) then
            do i=1,ni 
                  read(1,*) xlow(i), ylow(i), xhigh(i), yhigh(i)
            end do
      else
            print *, "ni or nj is too large"
      end if
      
! INSERT your code here

! Now read in the flow data from unit 2.
! You should read in the following variables sequentially:

      read(2,*) rgas, gamma
      read(2,*) pstagin, tstagin, alpha1, pdown
      read(2,*) cfl, smooth_fac_in
      read(2,*) nsteps, conlim_in
!       rgas, gamma
!       pstagin, tstagin, alpha1, pdown
!       cfl, smooth_fac_in
!       nsteps, conlim_in

! INSERT your code here

! Set some other variables that will be used throughout the
! calculation. Change alpha1 to radians. Scale the smoothing factor
! and convergence limits by cfl: changes over a timestep should be proportional
! to deltat (which is proportional to cfl).

      emax       = 1000000.0
      eavg       = emax
      cp         = rgas*gamma/(gamma-1.0)
      cv         = cp/gamma
      fga        = (gamma - 1.0)/gamma
      smooth_fac = smooth_fac_in*cfl
      conlim     = conlim_in*cfl
      alpha1     = alpha1*3.14159/180.0

! Close the flow and geom files

      close(1)
      close(2)

! Calculate the reference values which are used to check convergence

      ncells    =  ni * nj
      jmid      =  (1 + nj)/2
      roin      =  pstagin/rgas/tstagin
      ref_ro    =  (pstagin-pdown)/rgas/tstagin
      ref_t     =  tstagin*(pdown/pstagin)**fga
      ref_v     =  sqrt(2*cp*(tstagin-ref_t))
      ref_rovx  =  roin*ref_v
      ref_rovy  =  ref_rovx
      ref_roe   =  roin*cv*(tstagin-ref_t)

      end
      subroutine set_fluxes

      use common_block


      implicit none

! Local stuff
      integer ::  i, j

! This subroutine calculates the fluxes of mass momentum and energy
! across the faces of every cell.

! The "i" faces with i = 1 or i = ni are the upstream and downstream
! boundaries to the flow domain, while "j" faces with j = 1 or j = nj
! are solid boundaries. All fluxes are calculated assuming a linear
! variation in the flow properties between the cell corner nodes.

! First calculate the mass flux across each "i" face of the elements.
! Also calculate the total mass flow rate "flow(i)" across each "i" line.
! This will be used to check global continuity.

      do i=1,ni
        flow(i) = 0.0
        do j=1,nj-1
          fluxi_mass(i,j) = 0.5*( (rovx(i,j)+rovx(i,j+1))*dlix(i,j) +  &
                                  (rovy(i,j)+rovy(i,j+1))*dliy(i,j) )
          flow(i) = flow(i) + fluxi_mass(i,j)
        end do
      end do

! Now the mass flux across each "j" face.

      do i=1,ni-1
        do j=2,nj-1
          fluxj_mass(i,j) = 0.5*( (rovx(i,j)+rovx(i+1,j))*dljx(i,j) +  &
                                  (rovy(i,j)+rovy(i+1,j))*dljy(i,j) )
! INSERT your code here to calculate "fluxj_mass(i,j)"

        end do
      end do

! Set the mass fluxes through the j=1 and j=nj faces to zero as
! these are solid surfaces. It is not necessary to resolve the
! velocity parallel to the surfaces.

      do i=1,ni-1
        fluxj_mass(i,1) = 0.0
        fluxj_mass(i,nj)= 0.0
      end do

! Calculate the fluxes of x-momentum
! INSERT your code here to set "fluxi_xmom(i,j)"
      do i=1,ni
        do j=1,nj-1
          fluxi_xmom(i,j) =  0.5*fluxi_mass(i,j)*(vx(i,j)+vx(i,j+1)) &
                         + 0.5*(p(i,j)+p(i,j+1))*dlix(i,j)
          
          ! 0.5*( fluxi_mass(i,j)*(vx(i,j)+vx(i,j+1))  +  &
          !                       (p(i,j)+p(i,j+1))*dlix(i,j)   )
        end do
      end do

! INSERT your code here to set "fluxj_xmom(i,j)"
      do i=1,ni-1
        do j=1,nj
          fluxj_xmom(i,j) =  0.5*fluxj_mass(i,j)*(vx(i,j)+vx(i+1,j)) &
                         + 0.5*(p(i,j)+p(i+1,j))*dljx(i,j)
                    ! fluxj_xmom(i,j) = 0.5*( fluxj_mass(i,j)*(vx(i,j)+vx(i+1,j))  +  &
          !                       (p(i,j)+p(i+1,j))*dliy(i,j)   ) 
        end do
      end do

! Calculate the fluxes of y-momentum
! INSERT your code here to set "fluxi_ymom(i,j)"
      do i=1,ni
        do j=1,nj-1
          fluxi_ymom(i,j) =  0.5*fluxi_mass(i,j)*(vy(i,j)+vy(i,j+1)) &
                         + 0.5*(p(i,j)+p(i,j+1))*dliy(i,j)
          ! fluxi_ymom(i,j) = 0.5*( fluxi_mass(i,j)*(vy(i,j)+vy(i,j+1))  +  &
          !                       (p(i,j)+p(i,j+1))*dljx(i,j)   )
        end do
      end do

! INSERT your code here to set "fluxj_ymom(i,j)"
      do i=1,ni-1
        do j=1,nj
          fluxj_ymom(i,j) =  0.5*fluxj_mass(i,j)*(vy(i,j)+vy(i+1,j)) &
                         + 0.5*(p(i,j)+p(i+1,j))*dljy(i,j)
          ! fluxj_ymom(i,j) = 0.5*( fluxj_mass(i,j)*(vy(i,j)+vy(i+1,j))  +  &
          !                       (p(i,j)+p(i+1,j))*dljy(i,j)   )
        end do
      end do

! Calculate the fluxes of enthalpy
! INSERT your code here to set "fluxi_enth(i,j)"
      do i=1,ni
        do j=1,nj-1
          fluxi_enth(i,j) = 0.5*fluxi_mass(i,j)*(hstag(i,j)+hstag(i,j+1))
        end do
      end do
! INSERT your code here to set "fluxj_enth(i,j)"
      do i=1,ni-1
        do j=1,nj
          fluxj_enth(i,j) = 0.5*fluxj_mass(i,j)*(hstag(i,j)+hstag(i+1,j))
        end do
      end do
! Note that we could have set the flux of enthalpy to zero on
! j=1 and j=nj. This would save a bit of cpu time but the fluxes
! will be zero anyhow since the mass fluxes were set to zero.

      end
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
      subroutine set_timestep

      use common_block


      implicit none

! Local stuff
      real  ::  umax
      integer :: i,j
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
                  tstag_local = (hstag(i,j)+hstag(i+1,j)+hstag(i,j+1)+hstag(i+1,j+1))/(4*cp)
                  c_local = sqrt(gamma*rgas*(tstag_local))
                  v_local = (norm2([vx(i,j),vy(i,j)]) + &
                              norm2([vx(i+1,j),vy(i+1,j)]) + &
                              norm2([vx(i,j+1),vy(i,j+1)]) + &
                              norm2([vx(i+1,j+1),vy(i+1,j+1)]))/4
                  step(i,j) = cfl*dmin(i,j)/abs(c_local+v_local)
            enddo
         enddo
      end

subroutine smooth(prop, corr_prop)


      use common_block


      implicit none

! This subroutine smooths the variable "prop" (i.e. it adds the
! artificial viscosity) by taking (1-sf) x the calculated value of
! "prop" + sf x (the average of the surrounding values of "prop").
! where sf is the smoothing factor.
! It is modified if the deferred correction extension is adopted.

      real, dimension(i_max,j_max) ::  prop, corr_prop

! Local stuff
      real, dimension(i_max,j_max) ::  store
      real  ::  sf, sfm1, avg, avg1, avgnj, corrnew
      integer  ::  i, j, ip1, im1

! To avoid using already smoothed values for smoothing other values
! the smoothed values are initially stored in an array "store".
! Wall nodes do not have 4 neighbours, so the formula for smoothing them
! is different.
! Inlet and exit values are constrained by boundary conditions and extrapolation from
! previously smoothed nodes, so how they are smoothed has little material effect
! on the solution.  Remember that they do not have 4 neighbours.

      sf   = smooth_fac
      sfm1 = 1.0 - sf
      fcorr = 0.9

      do i=1,ni
         ip1 = i+1
         if( i==ni ) ip1 = ni
         im1 = i-1
         if( i==1  ) im1 = 1

         do j=2,nj-1
            avg = 0.25*(prop(ip1,j)+prop(im1,j)+prop(i,j-1)+prop(i,j+1))
            corrnew = fcorr*(prop(i,j)-avg)
            corr_prop(i,j) = 0.99*corr_prop(i,j) + 0.01*corrnew

! INSERT your code here
            store(i,j) = avg + corr_prop(i,j)
      !      store(i,j) = sfm1*prop(i,j) + sf*avg

         enddo

! On the surfaces j=1 and j=nj take the average as shown below.

         avg1  = (prop(im1,1)+prop(ip1,1)+2.*prop(i,2)-prop(i,3))/3.0
         avgnj = (prop(im1,nj) + prop(ip1,nj) + 2.*prop(i,nj-1)    &
                  - prop(i,nj-2))/3.0

! INSERT your code here to smooth the surface values
         store(i,1) = sfm1*prop(i,1) + sf*avg1
         store(i,nj) = sfm1*prop(i,nj) + sf*avgnj
      enddo

! Reset the smoothed value to "prop" before returning to the main program.

      do i=1,ni
        do j=1,nj
          prop(i,j) = store(i,j)
        end do
      end do

      end
      subroutine sum_fluxes(iflux,jflux,delprop,prop_inc, frkut_)

      use common_block


      implicit none

! This subroutine sums the fluxes for each element, calculates the changes
! in the variable "prop" appropriate to a cell (delprop) and distributes them to the
! four corners of the element (stored in prop_inc).

      real, dimension(i_max,j_max) ::  iflux, jflux, prop_inc, delprop,  &
                                       previous, store
     real :: frkut_

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
                  delprop(i,j) = (step(i,j)*frkut_/area(i,j))*(iflux(i,j) - iflux(i+1,j) + jflux(i,j) - jflux(i,j+1))
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
