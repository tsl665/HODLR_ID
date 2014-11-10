      subroutine LUBackSolveZ(N,NRHS,A,LDA,IPIV,B,LDB)

      integer N,NRHS,LDA,LDB
      complex *16 A(LDA,N),B(LDB,NRHS)
      integer IPIV(N)

      integer i,j,k
      complex *16 swap

c     First permute b using IPIV (invert P)
      do i = 1,N
cc          call prinf('IPIV(i) = *',IPIV(i),1)
cc          call prin2('B = *',B,N*NRHS*2)
            do j = 1, NRHS
                  swap = B(i,j)
                  B(i,j) = B(IPIV(i),j)
                  B(IPIV(i),j) = swap
            enddo !j
      enddo !i

c     Back Solve LY = P^-1 B

      do j = 1, NRHS
            do i = 2, N
                  do k = 1, i-1
                        B(i,j) = B(i,j) - A(i,k) * B(k,j)
                  enddo !k
            enddo !i
      enddo !j

c     Back Solve UX = Y

      do j = 1, NRHS
            do i = N, 1, -1
                  do k = N, i+1, -1
                        B(i,j) = B(i,j) - A(i,k) * B(k,j)
                  enddo !k
                  B(i,j) = B(i,j) / A(i,i)
            enddo !i
      enddo !j

      return
      end
