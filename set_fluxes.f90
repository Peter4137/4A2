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
          fluxi_xmom(i,j) = 0.5*( fluxi_mass(i,j)*(vx(i,j)+vx(i,j+1))  +  &
                                (p(i,j)+p(i,j+1))*dlix(i,j)   )
        end do
      end do

! INSERT your code here to set "fluxj_xmom(i,j)"
      do i=1,ni-1
        do j=1,nj
          fluxj_xmom(i,j) = 0.5*( fluxj_mass(i,j)*(vx(i,j)+vx(i+1,j))  +  &
                                (p(i,j)+p(i+1,j))*dliy(i,j)   ) 
        end do
      end do

! Calculate the fluxes of y-momentum
! INSERT your code here to set "fluxi_ymom(i,j)"
      do i=1,ni
        do j=1,nj-1
          fluxi_ymom(i,j) = 0.5*( fluxi_mass(i,j)*(vy(i,j)+vy(i,j+1))  +  &
                                (p(i,j)+p(i,j+1))*dljx(i,j)   )
        end do
      end do

! INSERT your code here to set "fluxj_ymom(i,j)"
      do i=1,ni-1
        do j=1,nj
          fluxj_ymom(i,j) = 0.5*( fluxj_mass(i,j)*(vy(i,j)+vy(i+1,j))  +  &
                                (p(i,j)+p(i+1,j))*dljy(i,j)   )
        end do
      end do

! Calculate the fluxes of enthalpy
! INSERT your code here to set "fluxi_enth(i,j)"
      do i=1,ni
        do j=1,nj-1
          fluxi_enth(i,j) = 0.5*( rovx(i,j)+rovx(i,j+1))*0.5*(hstag(i,j)+hstag(i,j+1)) + &
                            0.5*( rovy(i,j)+rovy(i,j+1))*0.5*(hstag(i,j)+hstag(i,j+1))
        end do
      end do
! INSERT your code here to set "fluxj_enth(i,j)"
      do i=1,ni-1
        do j=1,nj
          fluxj_enth(i,j) = 0.5*( rovx(i,j)+rovx(i+1,j))*0.5*(hstag(i,j)+hstag(i+1,j)) + &
                            0.5*( rovy(i,j)+rovy(i+1,j))*0.5*(hstag(i,j)+hstag(i+1,j))
        end do
      end do
! Note that we could have set the flux of enthalpy to zero on
! j=1 and j=nj. This would save a bit of cpu time but the fluxes
! will be zero anyhow since the mass fluxes were set to zero.

      end
