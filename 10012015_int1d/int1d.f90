module calcU
  double precision r1l, r2l
  double precision r12
  double precision gama, y, s
  contains
    function getU(x)
      ! x is lambda
      double precision getU, x
      double precision tmp
      tmp = y + 2.0d0 * s * x
      getU = erf(0.5d0 * tmp * (x * (1.0d0 - x))**(-0.5d0)) / tmp * gama
    end function getU
    
    function getU2(x)
      double precision getU2, x
      double precision tmp
      tmp = (1.0d0 - x) * r2l + x * r1l
      getU2 = 1.0d0 / tmp * erf(tmp * (0.5d0 / (2.0d0 * x * (1 - x)))**(0.5))
    end function getU2

    subroutine init(r1, r2, beta)
      double precision r1, r2, beta
      double precision r12, eief, mu, hbar
      r1l = r1
      r2l = r2
      r12 = (r2 + r1) / 2.0d0
      r12 = r2 - r1
      eief = 1.0d0
      mu = 0.5d0
      hbar = 1.0d0
      gama = eief * (2.0d0 * beta * mu)**0.5d0 / hbar
      y = (r1 + r2 - r12) * (2.d0 * hbar**2.d0 * beta / mu)**(-0.5d0)
      s = r12 * (2.d0 * hbar**2.d0 * beta / mu)**(-0.5d0)
    end subroutine
end module calcU

program main
  use calcU
  double precision integrate, f, a, b, eps
  double precision r
  integer maxLoop
  integer i
  external integrate
  external f
  a = 0.0d0 + 1.0d-6
  b = 1.0d0 - 1.0d-6
  eps = 1.0d-6
  maxLoop = int(1e5)

  if(0 < 1) then
  do i = 1, 30
    r = i * 1.d-1
    call init(r, r, 1.0d0)
    print *, r, r, integrate(getU2, a, b, eps, maxLoop)
  end do
  print *, "-----"
  do i = 1, 30
    r = i * 1.d-1
    call init(-r, r, dble(1.0))
    print *, -r, r, integrate(getU, a, b, eps, maxLoop)
  end do
  print *, "------"
  end if

  call init(1.0d0, 1.0d0, 1.0d0)
  print *, integrate(getU, a, b, eps, maxLoop)
end program main

function f(x)
  implicit none
  double precision f, x
  f =  x * x
end function f

function integrate(f, a, b, eps, maxLoop)
  ! Integrate f from a to b using rk45
  implicit none
  double precision integrate, f, a, b, EPS_DEFAULT
  double precision h, x, w, w1, w2, k1, k2, k3, k4, k5, k6, r, delta
  integer MAX_LOOP_DEFAULT, cnt
  double precision eps
  integer maxLoop
  
  w = 0.0d0 ! initial sum
  h = 1.0d-3 ! initial step size
  x = a ! integrate start from a
  
  cnt = 0
  do while (cnt < maxLoop)
    if (x >= b) then
      exit
    end if
    h = min(h, b - x)
    k1 = h * f(x)
    k2 = h * f(x + h * 0.25d0)
    k3 = h * f(x + h * 0.375d0)
    k4 = h * f(x + h * 12d0 / 13d0)
    k5 = h * f(x + h)
    k6 = h * f(x + h * 0.5d0)
    w1 = w + 25d0 * k1 / 216d0 + 1408d0 * k3 / 2565d0 + 2197d0 * k4 / 4104d0 - k5 / 5d0;
    w2 = w + 16d0 * k1 / 135d0 + 6656d0 * k3 / 12825d0 + 28561d0 * k4 / 56430d0 - 9d0 * k5 / 50d0 + 2d0 * k6 / 55d0
    r = abs(w1 - w2) / h

    ! my change to rk45 to avoid underflow r == 0
    if (r <= eps) then
      r = 1.d-60
    end if

    delta = 0.84d0 * (eps / r)**0.25d0
    if (r <= eps) then
      x = x + h
      w = w1
      cnt = cnt + 1
    end if
    h = delta * h
    !print *, "Step", cnt, x, w1, w2, delta, r, h
    !print *, h, k1, k2, k3, k4, k5, k6, w1, w2, r, delta, h, w
  end do

  integrate = w
end function integrate


