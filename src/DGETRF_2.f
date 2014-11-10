      subroutine ZGETRF_2(N,A,LDA,LU,IPIV)
c     Computes LU fatorization of square Matrix A without changing 
c     original matrix, and store the fatorization into new 
c     matrix LU.
      implicit none
      integer N, LDA
      complex*16 A(LDA,*), LU(LDA,N)
      integer IPIV(*)
      integer INFO, i, j


      do j = 1,N
            do i = 1,N
                  LU(i,j)=A(i,j)
            enddo
      enddo

      call ZGETRF(N,N,LU,LDA,IPIV,INFO)

      return
      end



