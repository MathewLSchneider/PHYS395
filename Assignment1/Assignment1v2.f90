! Assignment1.f90  -  starting point for Gauss-Jordan elimination
! compile with: gfortran -O3 -fdefault-real-8 -o GaussJ Assignment1.f90

program test
implicit none

integer  n, i, j
real F(100), x(100), T(100,100), C(100)
n = 100

do i = 1,n
	x(i) = -1.0 + (i-1.0)*2.0/(n-1)
end do
forall (i=1:n) F(i) = (1.0/(1.0 + 10.0*x(i)*x(i)))

do i=1,n
	do j=1,n
		T(i,j) = ChebyshevT(x(i),j-1)
	end do
end do

 call gaussj(n, T, F)
!do i = 1,n
!	write (*,*) x(i), C(i,:)
!end do
do i=1,n
	do j = 1,n
		C(i) = C(i) + F(j)*ChebyshevT(x(i),j)
	end do
end do
do i=1,n
	write (*,*) x(i), (1.0/(1.0 + 10.0*x(i)*x(i))), C(i)
end do

contains


! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
subroutine gaussj(n, A, B)
	integer n, i, j; real A(n,n), B(n), M(n,n+1)

!Make matrix M which is just [A|B]
	M(:,1:n) = A(:,:)
	M(:,n+1) = B(:)

	!Triangularize the A part of the matrix
	do j = 1,n-1
		do i = j+1,n
			M(i,:) = M(i,:) - (M(i,j)/M(j,j)) * M(j,:)
		end do
	end do

	!Diagonalize and normalize the A part of the matrix
	do j = n,2,-1
		do i = j-1,1,-1
			M(i,:) = M(i,:) - (M(i,j)/M(j,j)) * M(j,:)
		end do
		M(j,:) = M(j,:)/M(j,j)
	end do

	
	!The last column is now the solution
	B(:) = M(:,n+1)

end subroutine


! Chebyshev polynomial T_n(x)
elemental function ChebyshevT(x, n)
	real ChebyshevT, x; integer n
	intent(in) x, n
	
	ChebyshevT = cos(n*acos(x))
end function


end