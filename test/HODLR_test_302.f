      implicit none
      integer n, r, p, MAXLEN
c     parameter (MAXLEN = 1E7)
c     parameter (n = 35, p = 3, r = 1)
      real *8 done
      parameter (done = 1.0)
      complex *16 imag
      parameter (imag = (0,1))
      complex *16, allocatable :: A(:),AA(:),B(:),BB(:),C(:)
      real *8, allocatable :: error(:)
      real *8 ts,te, TIME(10)
c     time start, time end
      integer i,j,k

      print *, 'Enter n (order of HODLR Matrix):'
      read *, n

      print *, 'Enter p (rank of low-rank approximation):'
      read *, p


      print *, 'Enter r (Number of R.H.S.):'
      read *, r

      allocate(A(n*n))
      allocate(AA(n*n))
      allocate(B(n*r))
      allocate(BB(n*r))
      allocate(C(n*r))
      allocate(error(r))


      call prini(6,13)
      call generateHODLR(n,A)
      call copyMatrixZ(n,n,A,AA,n,n)
      
c     call generateRHS()
      do j = 1,r
            do i = 1,n
                  B(i+(j-1)*n) = cos(i*done)+imag*i
                  C(i+(j-1)*n) = B(i+(j-1)*n)
            enddo
      enddo
      call copyMatrixZ(n,r,B,BB,n,n)

c     call CPU_TIME(ts)

      call HODLR_ID(n,A,n,p,B,n,r,TIME)

c     call CPU_TIME(te)

      call prin2('Elapsed CPU time = *', TIME, 4)

      call bErrorZ(n,n,r,AA,BB,B,n,n,n,error)

      call prin2('Backward Error = *', error, r)



c     i = 1
c     j = 1
c     call passvalue(i)
c     call prinf('i - j =*',i-j,1)

      stop
      end

      subroutine passvalue(n)
      value n
      integer n
      n = n+1
      return
      end
      
      subroutine matrixMul(n,m,l,LDA,LDB,LDC,A,B,C)
c     matrix multiplication: C = A*B
c     A: m-by-l
c     B: l-by-n
c     C: m-by-n
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
c     Computing the backward error for AX = B
c     A: m-by-n (usually m = n)
c     B: m-by-r
c     X: n-by-r
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


c     call prin2('termwise error = *', C, m*r*2)

      do j = 1,r
            error(j) = 0
            do i = 1,m
                  error(j) = error(j) + C(i,j) * CONJG(C(i,j))
            enddo
            error(j) = sqrt(error(j))
      enddo

      return
      end



