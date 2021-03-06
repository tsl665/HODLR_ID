	implicit none
	integer nLEN
	parameter (nLEN = 10)
	integer n(nLEN),r,p
	parameter (r = 1, p = 5)
	real *8 time(nLEN)
c	real *8, allocatable :: error(:)
	integer i,j,k

	call prini(6,13)
	call prinf('r = *',r,1)
	call prinf('p = *',p,1)
	do i = 1, nLEN
		n(i) = 100*(2**(i-1))
		call prinf('n = *',n(i),1)
		call unitTest(n(i),p,r,time(i))
		call prin2('time = *', time(i), 1)
c		call prin2('error = *', error, r)
	enddo



	stop
	end


	

	





	subroutine unitTest(n,p,r,time)
	implicit none
	integer n, r, p, MAXLEN
c	parameter (MAXLEN = 1E7)
c	parameter (n = 35, p = 3, r = 1)
	complex *16, allocatable :: A(:),AA(:),B(:),BB(:),C(:)
	real *8, allocatable :: error(:)
	real *8 ts,te,time
c	time start, time end
	integer i,j,k

c	print *, 'Enter n (order of HODLR Matrix):'
c	read *, n

c	print *, 'Enter p (rank of low-rank approximation):'
c	read *, p


c	print *, 'Enter r (Number of R.H.S.):'
c	read *, r

	allocate(A(n*n))
	allocate(AA(n*n))
	allocate(B(n*r))
	allocate(BB(n*r))
	allocate(C(n*r))
	allocate(error(r))


c	call prini(6,13)
	call generateHODLR(n,A)
	call copyMatrixZ(n,n,A,AA,n,n)
	
c	call generateRHS()
	do j = 1,r
		do i = 1,n
			B(i+(j-1)*n) = 1
			C(i+(j-1)*n) = 1
		enddo
	enddo
	call copyMatrixZ(n,r,B,BB,n,n)

	call CPU_TIME(ts)

	call HODLR_ID(n,A,n,p,B,n,r)

	call CPU_TIME(te)
	time = te - ts

c	call prin2('Elapsed CPU time = *', te-ts, 1)

	call bErrorZ(n,n,r,AA,BB,B,n,n,n,error)

	call prin2('Backward Error = *', error, r)

	return
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




