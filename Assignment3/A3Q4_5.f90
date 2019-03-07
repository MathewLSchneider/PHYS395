! A3Q4_5.f90  -  Do nonlinear least squares fit using 
!										Levenberg-Marquardt
! compile with: gfortran -O3 -fdefault-real-8 A3Q4_5.f90 -llapack

program A3Q4_5
implicit none

integer i, n, k
real x, y, C(5), B(4), dfdc(5), Lambda(5,5), A(5,5), R(5), Chi(2)
real, parameter :: pi = 3.14152653589793238462643383279502884197
real, dimension(:,:), allocatable :: data, J, Fn
integer :: length, status, pivot(5), info

!Initialize a bunch of variables
length = 0; status = 0; B = 0; Lambda = 0; R = 0
forall (i=1:5) Lambda(i,i) = 1000.0

!Initial guess for coefficients
C = 0
C(2) = 2.4; C(5) = -5.3

!Get length of data
do while (status == 0)
	read (*,*,iostat=status) x, y
	if (status < 0) exit
	length = length + 1
end do

!Read in data
open (unit = 2, file='Assignment #2.dat', action='read', &
	status='old')
allocate(data(length,2)); 
data = 0; 
do i=1,length
	read (2,*) data(i,1), data(i,2)

end do
close(2)

allocate(J(length,5)); allocate(Fn(length,1))
Chi(1) = 0; Chi(2) = 1000

!Iterate until Chi squared does not change much
do while (abs(Chi(1) - Chi(2)) > 5.0)

	!Define the Jacobian to use in Levenberg-Marquardt method
	do k=1,5
		do i=1,length
			call Basis(data(i,1),B)
			dfdc = df(C, B)
			J(i,k) = dfdc(k)
		end do
	end do

	!Get function values
	do i=1,length
		call Basis(data(i,1),B)
		Fn(i,1) = f(C,B)
	end do

	!Calculate Chi squared
	Chi(1) = Chi(2)
	Chi(2) = sum(abs(data(:,2) - Fn(:,1))**2)
	
	!Calculate the matrices for the lin equation solve for L-M
	A = matmul(transpose(J),J) + Lambda
	R = matmul(transpose(J), data(:,2) - Fn(:,1))

	!Solve for new values of coefficients
	call dgesv(5, 1, A, 5, pivot, R, 5, info)
	C = C + R

end do

write (*,*) 'Fit done'

write (*,*) 'Chi squared is:', Chi(2)
write (*,*) 'Chi squared should be approximately', length
write (*,*) 'Coefficients are :', C
write (*,*) 'Where C1-C4 are in the exponent, C5 is the constant'

open (unit=3,file='FIT.txt',action='write', &
	status='new',position='rewind')
do i=1,length
	write (3,*) data(i,1), Fn(i,1)
end do
close(3)

contains

!Calculate value of basis functions at x
subroutine Basis(x, B)
	real x, B(4)

	B = [1.0, cos(2.0*pi*x), cos(4.0*pi*x), cos(6.0*pi*x)]

end subroutine

! function to fit
function f(C, B); 
	real f, C(5), B(4)
	f = exp(sum(C(1:4)*B)) + C(5)
end function

!Derivatives of function to fit with respect to coefficients
function df(C,B)
	real df(5), C(5), B(4)
	integer i

	forall (i=1:4) df(i) = B(i)*exp(sum(C(1:4)*B)) 
	df(5) = 1.0
end function




end program