	subroutine SMW_ID(A,LDA,work,List,Proj,B,LDB,r,pw)
c	MUST IMPROVE
c	In future version, should compute SMW inverse without
c	explicitly forming U and V.
	implicit none
	integer LDA, LDB, r, pw,List(*),work(*)
	complex *16 A(LDA,*),B(LDB,r),Proj(*)

	integer p,pwL,pwR,nL,nR,pp,pl,irow,icol,i,j,k
	integer IPVP(2*work(37)),INFO
	complex *16 U(work(pw+work(4)),2*work(37))
	complex *16 VT(2*work(37),work(pw+work(4)))
	complex *16 C(2*work(37),2*work(37))
	complex *16 y(2*work(37),r)
	complex *16 x(work(pw+work(4)),r)
	

	p = work(37)
	pwL = work(pw + work(7)) ! lc of pw
	pwR = work(pw + work(8)) ! rc of pw
	nL = work(pwL + work(4))
	nR = work(pwR + work(4))

c	Initialize U and VT by zero matrix

	do j = 1,2*p
		do i = 1,work(pw+work(4))
			U(i,j) = 0
		enddo !i
	enddo !j

	do j = 1,work(pw+work(4))
		do i = 1,2*p
			VT(i,j) = 0
		enddo !i
	enddo !j
	
c	Construct upper-left block of U (copies of submatrix of A)
	pl = work(pwL+work(10))
	do j = 1,p
		icol = List(pl+j-1)
		do i = 1,nL
			irow = work(pwL+work(3))+i-1
			U(i,j)=a(irow,icol)
		enddo !i
	enddo !j

c	Construct lower-right block of U
	pl = work(pwR+work(10))
	do j = 1,p
		icol = List(pl+j-1)
		do i = 1,nR
			irow = work(pwR+work(3))+i-1
			U(nL+i,p+j)=a(irow,icol)
		enddo !i
	enddo !j

	
c	Construct upper-right block of VT
	pp = work(pwL+work(11))
	do j=1,nR
		do i=1,p
			VT(i,j+nL) = Proj(pp+i+(j-1)*p-1)
		enddo
	enddo

c	Construct lower-left block of VT
	pp = work(pwR+work(11))
	do j=1,nL
		do i=1,p
			VT(i+p,j) = Proj(pp+i+(j-1)*p-1)
		enddo
	enddo


c	SMW Formula
	y = matmul(VT,b(1:work(pw+work(4)),1:r)) 
	!probably can't use matmul when b is submatrix
c	stop
	C = matmul(VT,U)
	do i = 1, 2*p
		C(i,i) = C(i,i)+1
	enddo


c	check here for LDB
	call ZGESV_2(2*p,r,C,2*p,y,2*p)

	x = matmul(U,y)

	do j = 1,r
		do i = 1,work(pw+work(4))
			b(i,j) = b(i,j)-x(i,j)
		enddo
	enddo

	return
	end
	

	subroutine ZGESV_2(N,NRHS,A,LDA,B,LDB)
	integer N, NRHS,LDA, LDB
	complex *16 A(LDA,*),B(LDB,*)
	integer IPIV(N), INFO

	complex *16 AA(LDA,N)
	integer i,j

	do i = 1,N
		do j = 1,N
			AA(i,j)=A(i,j)
		enddo
	enddo

	call ZGESV(N,NRHS,AA,LDA,IPIV,B,LDB,INFO)

	return
	end subroutine ZGESV_2
	
	subroutine ZGESV_3(N,NRHS,A,LDA,B,LDB,LU,IPIV)
	integer N, NRHS,LDA, LDB, IPIV(N)
	complex *16 A(LDA,*),B(LDB,*)
	integer INFO

	complex *16 LU(LDA,N)
	integer i,j

c	complex *16 LU2(N,N)
c	complex *16 E(N,N)
c	real *8 Err


	do i = 1,N
		do j = 1,N
			LU(i,j)=A(i,j)
c			LU2(i,j)=A(i,j)
c			E(i,j) = LU(i,j)-LU2(i,j)
		enddo
	enddo

c	do i = 1,N
c		call prin2('E_before(:,i) = *', E(1,i), N*2)
c	enddo
	call ZGESV(N,NRHS,LU,LDA,IPIV,B,LDB,INFO)
c	call ZGESV(N,NRHS,LU2,N,IPIV,B,LDB,INFO)

c	Err = 0

c	do i = 1,N
c		do j = 1,N
c			Err = Err + (LU(i,j)-LU2(i,j))**2
c		enddo
c	enddo

c	Err = sqrt(Err)
c	call prin2('err = *', err, 1)


	return
	end subroutine ZGESV_3

	
	subroutine checkLR(m,n,krank,A,LDA,List,Proj)
	integer m,n,krank,LDA
	complex *16 A(LDA,n), Proj(krank,n), B(m,krank), check(m,n)
	integer List(krank)

	integer i,j,k

	do j =1,krank
		do i = 1,m
			B(i,j) = A(i,List(j))
		enddo
	enddo
	
	check = matmul(B,Proj)

	do j = 1,n
		do i = 1,m
			check(i,j)=check(i,j)-A(i,j)
		enddo
	enddo

	call prin2('checkLR = *',check,n*m*2)

	

	return
	end

