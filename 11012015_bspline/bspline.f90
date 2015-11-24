module bspline

  use bspline_module
  use, intrinsic :: iso_fortran_env, only: wp => real64
  
  implicit none
  
  private

  integer :: nx, ny
  integer :: kx, ky

  real(wp) :: val
  integer :: i, j, k, idx, idy, iflag
  integer :: inbvx, inbvy
  integer :: iloy
  real(wp), allocatable :: tx(:), ty(:), bcoe(:, :)

  public :: bspline_init, bspline_eval

  contains
  
  ! x: x coordinates
  ! nx: number of x coordinates
  ! y: y coordinates
  ! ny: number of y coordinates
  ! data: functions values at each (x, y) grid point
  ! kx: x direction spline order
  ! ky: y direction spline order
  ! return: true or false, for success / fail, respectively
  subroutine bspline_init(x, nx_in, y, ny_in, fn_val, kx_in, ky_in)
    implicit none
    real(wp) :: x(:), y(:), fn_val(:, :)
    integer :: nx_in, ny_in, kx_in, ky_in
    ! logical :: bspline_init
    nx = nx_in
    ny = ny_in
    kx = kx_in
    ky = ky_in
    allocate(tx(nx + kx))
    allocate(ty(ny + ky))
    allocate(bcoe(nx, ny))
    idx = 0
    idy = 0
    inbvx = 1
    inbvy = 1
    iloy = 1
    call db2ink(x, nx, y, ny, fn_val, kx, ky, tx, ty, bcoe, iflag)
    ! bspline_init = iflag
  end subroutine bspline_init
  
  ! x, y: coordinate of point to evaluate spline interpolated value
  function bspline_eval(x, y)
    implicit none
    real(wp) :: x, y, bspline_eval
    call db2val(x, y, idx, idy, tx, ty, nx, ny, kx, ky, bcoe, val, iflag, inbvx, inbvy, iloy)
    bspline_eval = val
  end function bspline_eval

end module bspline
