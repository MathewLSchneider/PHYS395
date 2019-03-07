! A3Q1.f90  -  find a 3 roots by Newton's method
! compile with: gfortran -O3 -fdefault-real-8 A3Q1.f90

program A3Q1
implicit none

real x(3)
integer i,n
real, parameter :: phi = (1.0 + sqrt(5.0))/2.0
real, parameter :: epsilon = 1.0e-8
real a, b, c, d
! initial guess
x = [1.0, 0.0, -1.0]
a = -1.5; b = -0.5; c = b - (b-a)/phi
d = a + (b-a)/phi
! Newton's iteration
do n = 1,3
	do i = 1,8
		x(n) = x(n) - f(x(n))/df(x(n))
		!write (*,*) x, f(x)
	end do
end do

do n = 1,3
	write (*,*) x(n), f(x(n))
end do
write (*,*) '~~~~~~~~~~~~~~~~~~~~~'
call bracket(a,b,c,d)
write (*,*) (c+d)/2.0

a = 0.5; b = 1.0; c = b - (b-a)/phi
d = a + (b-a)/phi

call bracket(a,b,c,d)
write (*,*) (c+d)/2.0

contains

! function to find a root of...
pure function f(x); intent(in) x
	real f, x
	f = x**3 - x + 0.25
end function

! derivative of a function to find a root of
pure function df(x); intent(in) x
	real df, x
	df = 3*x**2 - 1
end function

! Bracketing subroutine
subroutine bracket(a, b ,c , d)
	real a, b, c, d

	do while (abs(c-d) > epsilon)
		if (g(c) < g(d)) then; b = d; end if
		if (g(d) < g(c)) then; a = c; end if
		c = b - (b-a)/phi
		d = a + (b-a)/phi
	end do

end subroutine

!function to bracket
pure function g(x); intent(in) x
	real g, x
	g = (x**2 - 1)**2 + x
end function

end program