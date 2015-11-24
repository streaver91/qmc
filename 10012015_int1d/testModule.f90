module testModule
  double precision a

  contains
    function f(x)
      double precision f, x
      f = x * a
    end function f

    subroutine init()
      a = 1.0
    end subroutine init 
end module testModule

function ff(f)
  double precision ff, f
  !print *, dble(3.0)
  ff = f(dble(3.0))
end function ff

program test
  use testModule
  double precision x, y
  double precision ff
  external ff
  x = 2.0
  call init()
  print *, f(x)
  print *, ff(f)
end program test
