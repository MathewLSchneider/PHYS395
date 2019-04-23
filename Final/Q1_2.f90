! A3Q1_3.f90  -  find a 3 roots by Newton's method,
!						 then 2 minima by bracketing
! compile with: gfortran -O3 -fdefault-real-8 Q1_2.f90

program Q1_2
implicit none

real(8) x(3)
integer i,n
real, parameter :: phi = (1.0 + sqrt(5.0))/2.0
real, parameter :: epsilon = 1.0e-12
real(16) a, b, c, d
! initial guess
x = [1.0, -2.0, -4.0]

!Initial bracket edges and middle point values
!c,d chosen according to golden ratio
a = -3.5; b = -2.5; c = b - (b-a)/phi
d = a + (b-a)/phi


! Newton's iteration
do n = 1,3
	do i = 1,10
		x(n) = x(n) - f(x(n))/df(x(n))
		!write (*,*) x, f(x)
	end do
end do
write (*,*) 'Roots of cos(x) = x/5 are at x = '
do n = 1,3
	write (*,*) x(n)
end do
write (*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

write (*,*) 'Minima of x^4 + 3x^3 - 4x^2 - 3x + 4 are at:'
call bracket(a,b,c,d, 1.0)

write (*,*) g((c+d)/2.0), 'at x =', (c+d)/2.0
write (*,*) ''
a = 0.0; b = 1.5; c = b - (b-a)/phi
d = a + (b-a)/phi

write (*,*) 'And'
write (*,*) ''
call bracket(a,b,c,d, 1.0)

write (*,*) g((c+d)/2.0), 'at x =', (c+d)/2.0
write (*,*) ''
a = -1.0; b = 0.0; c = b - (b-a)/phi
d = a + (b-a)/phi

write (*,*) 'Maxima is'
write (*,*) ''
call bracket(a,b,c,d, -1.0)

write (*,*) g((c+d)/2.0), 'at x =', (c+d)/2.0

contains

! function to find a root of...
pure function f(x); intent(in) x
	real f, x
	f = cos(x) - x/5.0
end function

! derivative of a function to find a root of
pure function df(x); intent(in) x
	real df, x
	df = -sin(x) - 1.0/5.0
end function

! Bracketing subroutine
subroutine bracket(a, b ,c , d, sign)
	real(16) a, b, c, d
	real(8) sign

!Shrink the interval while it is larger than some tolerance
	do while (abs(c-d) > epsilon)

		!Decrease inteval size based on function value of middle points
		if (sign > 0.0) then
			if (g(c) < g(d)) then; b = d; end if
			if (g(d) < g(c)) then; a = c; end if
		end if
		if (sign < 0.0) then
			if (g(c) > g(d)) then; b = d; end if
			if (g(d) > g(c)) then; a = c; end if
		end if
		!Reset middle points such that they keep the golden ratio proportion
		c = b - (b-a)/phi
		d = a + (b-a)/phi
	end do

end subroutine

!function to bracket
pure function g(x); intent(in) x
	real(16) g, x
	g = x**4 + 3.0*x**3 -4.0*x*x -3.0*x +4.0
end function

end program