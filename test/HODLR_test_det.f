	implicit none
	integer n, r, p, MAXLEN
c	parameter (MAXLEN = 1E7)
c	parameter (n = 35, p = 3, r = 1)
	complex *16, allocatable :: A(:),AA(:),C(:)
	real *8 error
	real *8 ts,te, TIME(10),de_true,de
c	time start, time end
	integer i,j,k
	integer, allocatable :: IPIV(:)

	print *, 'Enter n (order of HODLR Matrix):'
	read *, n

	print *, 'Enter p (rank of low-rank approximation):'
	read *, p

	allocate(A(n*n))
	allocate(AA(n*n))
	allocate(C(n*r))
	allocate(IPIV(n))

	de_true = 1.0
	call prini(6,13)
	call generateHODLR(n,A)
c	call prin2('A = *',A,n*n*2)

c	do i = 1,n
c		do j = 1,n
c			if (i >=j) then
c				A((i-1)*n+j) = n-i
c			else
c				A((i-1)*n+j) = j+i
c			endif
c		enddo
c	enddo
	call ZGETRF_2(n,A,n,AA,IPIV)
c	call prinf('IPIV = *',IPIV,n)
	do i = 1,n
		de_true = de_true * AA((i-1)*n+i)
		if (IPIV(i)/=i) then
			de_true = - de_true
		endif
	enddo

	
	call prin2('det_true = *',de_true, 1)
c	call CPU_TIME(ts)

	call HODLR_ID_det(de,n,A,n,p,TIME)

c	call CPU_TIME(te)

c	call prin2('Elapsed CPU time = *', TIME, 3)

	error = (de_true - de)/de_true

	call prin2('Backward Error = *', error, 1)
c	call prin2('det_true = *',de_true, 1)
	call prin2('det = *',de,1)

	deallocate(A)
	deallocate(AA)
	deallocate(C)
	deallocate(IPIV)


	stop
	end

	subroutine matrixMul(n,m,l,LDA,LDB,LDC,A,B,C)
c	matrix multiplication: C = A*B
c	A: m-by-l
c	B: l-by-n
c	C: m-by-n
	implicit none
	integer m,n,l,LDA,LDB,LDC
	complex *16 A(LDA,l), B(LDB,n), C(LDC,n)
	integer i,j,k

	do j = 1,n
		do i = 1,m
			C(i,j) = 0
			do k = 1,l
				C(i,j) = C(i,j)+A(i,k)*B(k,j)
			enddo
		enddo
	enddo

	return
	end

	subroutine copyMatrixZ(m,n,A,B,LDA,LDB)
	integer m,n,LDA,LDB
	complex *16 A(LDA,n), B(LDB,n)

	integer i,j

	do i = 1,m
		do j = 1,n
			B(i,j) = A(i,j)
		enddo
	enddo

	return
	end

	subroutine bErrorZ(m,n,r,A,B,X,LDA,LDB,LDX,error)
c	Computing the backward error for AX = B
c	A: m-by-n (usually m = n)
c	B: m-by-r
c	X: n-by-r
	integer m,n,r,LDA,LDB,LDX
	complex *16 A(LDA,n),B(LDB,r),X(LDX,r)
	complex *16 C(m,r), z
	real *8 error(r)

	integer i,j,k

	C = matmul(A(1:m,1:n),X(1:n,1:r))
	do j = 1,r
		do i = 1,m
			C(i,j) = C(i,j) - B(i,j)
		enddo
	enddo


c	call prin2('termwise error = *', C, m*r*2)

	do j = 1,r
		error(j) = 0
		do i = 1,m
			error(j) = error(j) + C(i,j) * CONJG(C(i,j))
		enddo
		error(j) = sqrt(error(j))
	enddo

	return
	end



