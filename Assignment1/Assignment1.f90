! Assignment1.f90  -  starting point for Gauss-Jordan elimination
! compile with: gfortran -O3 -fdefault-real-8 -o GaussJ Assignment1.f90

program test
implicit none

real A(3,3), B(3)

! matrix to invert
A(1,:) = [1.0, 2.0, 3.0]
A(2,:) = [3.0, 2.0, 1.0]
A(3,:) = [0.0, 1.0, 0.0]

B = [1.0, 2.0, 3.0]

write (*,*) A(1,:)
write (*,*) A(2,:)
write (*,*) A(3,:)

write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

 call gaussj(3, A, B)

!write (*,*) B

contains

! solve A.x = B using Gauss-Jordan elimination
! A gets destroyed, answer is returned in B
subroutine gaussj(n, A, B)
	integer n, i, j; real A(n,n), B(n), M(n,n+1)

!Make matrix M which is just [A|B]
	M(:,1:n) = A(:,:)
	M(:,n+1) = B(:)

	write (*,*) M(1,:)
	write (*,*) M(2,:)
	write (*,*) M(3,:)

	write (*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

	!Triangularize the A part of the matrix
	do j = 1,n-1
		do i = j+1,n
			M(i,:) = M(i,:) - (M(i,j)/M(j,j)) * M(j,:)
		end do
	end do

	write (*,*) M(1,:)
	write (*,*) M(2,:)
	write (*,*) M(3,:)
	
	write (*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	
	!Diagonalize and normalize the A part of the matrix
	do j = n,2,-1
		do i = j-1,1,-1
			M(i,:) = M(i,:) - (M(i,j)/M(j,j)) * M(j,:)
		end do
		M(j,:) = M(j,:)/M(j,j)
	end do

	write (*,*) M(1,:)
	write (*,*) M(2,:)
	write (*,*) M(3,:)
	
	write (*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
	
	!The last column is now the solution
	B(:) = M(:,n+1)

	write (*,*) B
end subroutine

end