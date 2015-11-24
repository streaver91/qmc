program mytest

    use bspline_module, only: db2ink, db2val 
    use, intrinsic :: iso_fortran_env, only: wp => real64

    implicit none

    integer, parameter :: nx = 6
    integer, parameter :: ny = 6
    
    integer, parameter :: kx = 3
    integer, parameter :: ky = 3
    
    real(wp) :: x(nx), y(ny)
    real(wp) :: tx(nx + kx), ty(ny + ky)
    real(wp) :: fcn_2d(nx, ny), bcoe(nx, ny)
    
    real(wp) :: tol
    real(wp), dimension(6) :: val, tru, err, errmax
    
    logical :: fail
    integer :: i, j, k, l, m, n, idx, idy, iflag
    integer :: inbvx, inbvy
    integer :: iloy
    
    fail = .false.
    tol = 1.0e-14_wp
    idx = 0
    idy = 0
    
    inbvx = 1
    inbvy = 1
    iloy = 1
    
    ! set grid
    do i = 1, nx
      x(i) = dble(i - 1) / dble(nx - 1)
    end do
    do j = 1, ny
      y(j) = dble(j - 1) / dble(ny - 1)
    end do
    
    
    
    ! calculate value at grid points
    do i = 1, nx
      do j = 1, ny
        fcn_2d(i, j) = f2(x(i), y(j))
      end do
    end do
    
    ! interpolate
    iflag = 0
    call db2ink(x, nx, y, ny, fcn_2d, kx, ky, tx, ty, bcoe, iflag)
    
    ! evaluate error
    errmax = 0.0_wp
    do i = 1, nx
      do j = 1, ny
        call db2val(x(i), y(j), idx, idy, tx, ty, nx, ny, kx, ky, bcoe, val(2), iflag, inbvx, inbvy, iloy)
        tru(2) = f2(x(i), y(j))
        err(2) = abs(tru(2) - val(2))
        errmax(2) = max(err(2), errmax(2))
        write (*, *) i, j, x(i), y(j), val(2), tru(2), err(2)
      end do
    end do
    
    i = 2
    write(*, *) i, 'D: max error: ', errmax(i)
    
    contains
    
    real(wp) function f2(x, y)
      implicit none
      real(wp) :: x, y
      f2 = 0.5_wp * (x * exp(-x) + y * exp(-y))
    end function f2

end program mytest
