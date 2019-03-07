! A3Q1_3.f90  -  find a 3 roots by Newton's method,
!						 then 2 minima by bracketing
! compile with: gfortran -O3 -fdefault-real-8 A3Q1_3.f90

program A3Q1_3
implicit none

real x(3)
integer i,n
real, parameter :: phi = (1.0 + sqrt(5.0))/2.0
real, parameter :: epsilon = 1.0e-8
real a, b, c, d
! initial guess
x = [1.0, 0.0, -1.0]

!Initial bracket edges and middle point values
!c,d chosen according to golden ratio
a = -1.5; b = -0.5; c = b - (b-a)/phi
d = a + (b-a)/phi


! Newton's iteration
do n = 1,3
	do i = 1,8
		x(n) = x(n) - f(x(n))/df(x(n))
		!write (*,*) x, f(x)
	end do
end do
write (*,*) 'Roots of x^3 - x + 0.25 are at x = '
do n = 1,3
	write (*,*) x(n)
end do
write (*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

write (*,*) 'Minima of (x^2 - 1)^2 + x are at:'
call bracket(a,b,c,d)

write (*,*) g((c+d)/2.0), 'at x =', (c+d)/2.0
write (*,*) ''
a = 0.5; b = 1.0; c = b - (b-a)/phi
d = a + (b-a)/phi

write (*,*) 'And'
write (*,*) ''
call bracket(a,b,c,d)

write (*,*) g((c+d)/2.0), 'at x =', (c+d)/2.0

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

!Shrink the interval while it is larger than some tolerance
	do while (abs(c-d) > epsilon)

		!Decrease inteval size based on function value of middle points
		if (g(c) < g(d)) then; b = d; end if
		if (g(d) < g(c)) then; a = c; end if
		
		!Reset middle points such that they keep the golden ratio proportion
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