

program Assignment2_third
implicit none

! number of basis functions
integer, parameter :: n = 8

real x, y, A(n,n), B(n), U(n), Chi, c
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

call svdsolve(n, A, U, c, 1.0e-6)

write (*,*) 'The condition number is:'
write (*,*) c

! output best fit parameters
write (*,*) 'The coefficients are:'
write (*,*) U

!Write coefficients to file
open (unit=1,file='Fit.txt',action='write', &
	status='old',position='rewind')
write (1,*) s, U
close(1)

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
subroutine svdsolve(n, A, B, c, epsilon)
	integer n, status; real A(n,n), B(n), epsilon
	
	! SVD matrices
	real U(n,n), V(n,n), S(n), W(6*n), c
	
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
	
	c = maxval(S)/minval(S)
	
	! compute the solution using pseudo-inverse
	B = matmul(transpose(U),B)
	where (S > epsilon*S(1)); B = B/S; elsewhere; B = 0.0; end where
	B = matmul(transpose(V),B)
end subroutine

end program
