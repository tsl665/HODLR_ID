      subroutine generateHODLR(n,A)
      integer n
c     Dimension of HOLDER matrix
      complex *16 A(n,n)
c     Output matrix
      integer i,j
      
      real *8 x(n)
      
      call srand(1)
      call prini(6,13)

      do i = 1,n
c           x(i) = rand(0)*2-1
            x(i) = -1 + 2.0*(i-1)/(n-1)
      enddo


      do j = 1,n
            do i = 1,n
                  A(i,j) = exp(-(x(i)-x(j))*(x(i)-x(j)))
                  if (i == j) then
                        A(i,j)=A(i,j)+1
                  endif
            enddo
      enddo
      


      return
      end subroutine generateHODLR

      subroutine generateRHS(n,r,B)
      implicit none
      integer n,r
      complex *16 B(n,r)
      integer i,j
      real *8 done
      complex *16 ima
      parameter (done = 1.0,ima = (0,1))

      do j = 1,r
            do i = 1,n
                  B(i,j) = cos(i*done) + ima * sin(j*done)
            enddo
      enddo

      return
      end subroutine generateRHS

