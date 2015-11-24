program test
  implicit none
  
  double precision  f, x, a,g
  external f,g
  double precision, allocatable :: args(:)
  


  a = 0.0
  !allocate(args(2))
  !args(1:2) = (0.0, 1.0)
  
  interface
  function f(x, args)
  implicit none
  double precision :: f, x
  double precision, allocatable :: args(:)
  end function f
  end interface

  !print *, g(args)
  !print *, (g(f, args))
  print *, f(a, args)
end program test


function f(x, args)
  implicit none
  double precision :: f, x
  double precision, allocatable :: args(:)
  f = x * (args(1) + args(2))
  
end function f

!function g(args)
!  double precision g
!  double precision :: args(:)
!  g = f((0.0), args)+ f(dble(1.0), args)+ f(dble(2.0), args)
!end function g
function g(f, args)
  implicit none
  double precision g, f
  double precision dummy
  double precision, allocatable :: args(:)
  dummy = 1.0
  g = (f(dummy, args))
end function g

