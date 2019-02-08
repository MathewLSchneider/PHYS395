! general linear least squares fit
! compile with gfortran -O3 -fdefault-real-8 leastsq.f90 -llapack

program leastsq
implicit none

! number of basis functions
integer, parameter :: n = 4

real x, y, A(n,n), B(n), U(n), Chi
real, dimension(:), allocatable :: fit, xp

integer i, j
integer :: s = 0
integer :: status = 0
integer :: stat = 0

! initialize accumulators
A = 0.0
U = 0.0
! read data from standard input
do while (status == 0)
	read (*,*,iostat=status) x, y
	if (status < 0) exit
	! evaluate basis functions at point x
	call evalb(x, B)
	!write (*,*) B
	s = s + 1
	! accumulate least square problem matrices
	forall (i=1:n,j=1:n) A(i,j) = A(i,j) + B(i)*B(j); U = U + y*B
end do
Chi = 0
allocate(fit(s))
fit = 0
allocate(xp(s))
xp = 0
! output accumulated matrices
!write (*,*) A
!write (*,*) U
! solve for best fit parameters
!call lsolve(n, A, U)
write (*,*) 'The condition number is:'
write (*,*) maxval(abs(A))/minval(abs(A))
call svdsolve(n, A, U, 1.0e-6)
write (*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
! output best fit parameters
s = 0
do while (stat == 0)
	read (*,*,iostat=stat) x, y
	if (stat < 0) exit
	
	! evaluate basis functions at point x
	!write (*,*) B
	! accumulate least square problem matrices
	s = s + 1
	xp(s) = x
	do i=1,n 
		fit(s) = fit(s) + U(i)*cos(2.0*(i-1)*355.0*x/113.0)
	end do
	Chi = Chi+(fit(s)-y)**2
end do
write (*,*) 'The coefficients are:'
write (*,*) U
write (*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
write (*,*) 'Chi Squared value is:'
write (*,*) Chi
if (Chi < 3) Then
	write (*,*) 'Fit looks good.'
end if
	
open (unit=1,file='Fit.txt',action='write', &
	status='new',position='rewind')
write (1,*) s, U
contains

! basis functions we are fitting
subroutine evalb(x, B)
	real x, B(n)
	
	! degenerate to illustrate SVD
	
	forall (i=0:n-1) B(i+1) = cos(2.0*i*x*355.0/113.0)
end subroutine

! solve A.x = B using LAPACK xGESV
! A gets destroyed, answer is returned in B
subroutine lsolve(n, A, B)
	integer n, pivot(n), status; real A(n,n), B(n)
	
	! initialize status to all clear
	status = 0
	
	! call appropriate version of xGESV
	select case (kind(A))
		case(4); call sgesv(n, 1, A, n, pivot, B, n, status)
		case(8); call dgesv(n, 1, A, n, pivot, B, n, status)
		case default; call abort
	end select
	
	! abort at first sign of trouble
	if (status /= 0) stop "singular matrix in lsolve()" 
end subroutine

! solve A.x = B using LAPACK xGESVD
! A gets destroyed, answer is returned in B
subroutine svdsolve(n, A, B, epsilon)
	integer n, status; real A(n,n), B(n), epsilon
	
	! SVD matrices
	real U(n,n), V(n,n), S(n), W(6*n)
	
	! initialize status to all clear
	status = 0
	
	! call appropriate version of xGESV
	select case (kind(A))
		case(4); call sgesvd('A', 'A', n, n, A, n, S, U, n, V, n, W, 6*n, status)
		case(8); call dgesvd('A', 'A', n, n, A, n, S, U, n, V, n, W, 6*n, status)
		case default; call abort
	end select
	
	! abort at first sign of trouble
	if (status /= 0) stop "singular matrix in svdsolve()" 
	
	! compute the solution using pseudo-inverse
	B = matmul(transpose(U),B)
	where (S > epsilon*S(1)); B = B/S; elsewhere; B = 0.0; end where
	B = matmul(transpose(V),B)
end subroutine

end program
