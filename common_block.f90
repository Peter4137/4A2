
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
