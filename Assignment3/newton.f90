! A3Q1.f90  -  find a 3 roots by Newton's method
! compile with: gfortran -O3 -fdefault-real-8 A3Q1.f90

program A3Q1
implicit none

real x(3)
integer i,n

! initial guess
x = [1.0, 0.0, -1.0]

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

end program