program bspline_test
  
  use bspline_oo_module
  ! use bspline
  use, intrinsic :: iso_fortran_env, only: wp => real64

  implicit none

  ! number of points in each direction
  integer, parameter :: nx = 100
  integer, parameter :: ny = 100
  
  ! order of spline curve in each direction
  integer, parameter :: kx = 3
  integer, parameter :: ky = 3
  
  real(wp) :: x(nx), y(ny) ! grid points
  real(wp) :: fn_val(nx, ny) ! function value for interpolation
  real(wp) :: val_eval, val_true, err, err_max
  integer :: i, j
  integer :: idx, idy
  integer :: iflag
  type(bspline_2d) :: s2

  idx = 0
  idy = 0

  do i = 1, nx
    x(i) = dble(i - 1) / dble(nx - 1)
  end do

  do j = 1, ny
    y(j) = dble(j - 1) / dble(ny - 1)
  end do

  do i = 1, nx
    do j = 1, ny
      fn_val(i, j) = fn(x(i), y(j))
    end do
  end do

  write (*, *) 'Initialize'


  call s2%initialize(x, y, fn_val, kx, ky, iflag)
  ! call bspline_init(x, nx, y, ny, fn_val, kx, ky)
  
  write (*, *) 'Evaluate'

  err_max = 0.0_wp

  do i = 1, nx
    do j = 1, ny
      call s2%evaluate(x(i), y(j), idx, idy, val_eval, iflag)
      val_true = fn_val(i, j)
      err = val_eval - val_true
      err_max = max(err_max, abs(err))
      write (*, *) i, j, x(i), y(j), val_eval, val_true, err
    end do
  end do

  write (*, *), 'Max Error: ', err_max

  contains

  real(wp) function fn(x, y)
    implicit none
    real(wp) :: x, y
    fn = exp(-x) * exp(-y)
  end function fn

end program bspline_test
