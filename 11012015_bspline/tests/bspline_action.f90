program bspline_action
  
  use bspline, only: bspline_init, bspline_eval
  ! use bspline
  use, intrinsic :: iso_fortran_env, only: wp => real64

  implicit none

  ! number of points in each direction
  integer, parameter :: nx = 101
  integer, parameter :: ny = 21
  
  ! order of spline curve in each direction
  integer, parameter :: kx = 3
  integer, parameter :: ky = 3
  
  integer, parameter :: line_skip = 5 ! # of lines skipped at the beginning
  integer, parameter :: interpolate_points = 10 ! # of points interpolated between points

  real(wp) :: x(nx), y(ny) ! grid points
  real(wp) :: fn_val(nx, ny) ! function value for interpolation
  real(wp) :: val_eval, val_true, err, err_max
  integer :: i, j
  character(100) :: line
  real(wp) :: q, sq, act

  ! skip first few lines that are not data
  do i = 1, line_skip
    read (*, *) line
    ! write (*, *) line
  end do

  do i = 1, nx
    do j = 1, ny
      ! save to tmp var
      read (*, *) q, sq, act
      fn_val(i, j) = act
      ! obtain y (s over q) grid
      if (i == 1) then
        y(j) = sq
      end if
    end do
    ! obtain q grid
    x(i) = q
  end do

  call bspline_init(x, nx, y, ny, fn_val, kx, ky)
  
  ! write (*, *) 'Evaluate'

  err_max = 0.0_wp
   
  do i = 1, nx * interpolate_points
    do j = 1, ny * interpolate_points
      q = (i - 1) * x(nx) / (nx * interpolate_points)
      sq = (j - 1) * y(ny) / (ny * interpolate_points)
      act = bspline_eval(q, sq)
      write (*, *) q, sq, act
      ! write (*, *) i, j, x(i), y(j), val_eval, val_true, err
    end do
  end do

  ! write (*, *), 'Max Error: ', err_max

  contains

  real(wp) function fn(x, y)
    implicit none
    real(wp) :: x, y
    fn = exp(-x) * exp(-y)
  end function fn

end program bspline_action
