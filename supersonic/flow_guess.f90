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

      machlim = 10.0 !supersonic limit
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
!Supersonic initial guess
      do i=1,ni
            v_guess(i) = vdown
            ro_guess(i) = rodown
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
