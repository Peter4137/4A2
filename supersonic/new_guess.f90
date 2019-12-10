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

    do j=1,nj
    do i=1,ni-1
        dx  = x(i+1,jmid) - x(i,jmid)
        dy  = y(i+1,jmid) - y(i,jmid)
        ds  = sqrt(dx*dx + dy*dy)
        vlocal    = vinlet  + (vdown-vinlet)*float(i-1)/float(ni-1)
        rolocal   = roin    + (rodown-roin)*float(i-1)/float(ni-1)
        tlocal    = tinlet  + (tdown-tinlet)*float(i-1)/float(ni-1)
        
        vx(i,j)      = vlocal*dy/ds
        vy(i,j)      = vlocal*dx/ds
        rovx(i,j) = rolocal*xvel
        rovy(i,j) = rolocal*yvel
        ro(i,j)   = rolocal
        roe(i,j)  = rolocal*(cv*tlocal + 0.5*vlocal*vlocal)
    end do
    vx(i,j) = vdown*dx/ds
    vy(i,j) = vdown*dy/ds

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