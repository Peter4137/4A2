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
      call flow_guess
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
        nrkut_max = 5
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

            call smooth(ro)
            call smooth(rovx)
            call smooth(rovy)
            call smooth(roe)
        end do 
! Check convergence and write out summary every 5 steps

        if(mod(nstep,5)==0) then
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
