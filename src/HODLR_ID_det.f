c     V302: 
c     1. Using standard back solve to avoid unnecessary ZGESV
c     2. Add new return variable TIME
c     3. Put small subroutines into a seperate file (SMW_SGESV2.f)

      subroutine HODLR_ID_det(de,n,A,LDA,p,TIME)
c     n: size of matrix A
c     A: HODLR matrix to solve
c     LDA: Leading Dimension of A
c     p: specific rank
c     b: r.h.s
c     LDB: Leading Dimension of b
c     r: number of r.h.s
c     TIME: an real array that records Elapsed CPU time for each
c           step.
c     TIME(1): total time
c     TIME(2): construct tree
c     TIME(3): dense invert
c     TIME(4): low-rank invert
      implicit none
      integer WLEN
      parameter (WLEN = 1E7) ! Length of work
      integer n, LDA, p, r
      real *8 TIME(10),ts,te,de
      complex *16 A(LDA,n)
      
      integer work(WLEN)
c     work: the object stores the tree structure of HODLR matrix.
c     work(1): starting index of tree structure. The first 1 to 
c           work(1)-1 memories are for constants. Here work(1) = 51.
c     
c     work(2) - work(20) is constants about attribution of each
c     single node. Let "tag" be the information we want. "pw"
c     is pointer for current node. Then we access that information
c     by work(pw+tag). "tag" is one of the following:
c
c     work(2): [lvl] Level of current node
c     work(3): [mstart] The starting column(row) of current node
c           in A
c     work(4): [msize] Size of current node
c     work(5): [odsci] Off-Diagonal Starting Column Index
c           Notice that the [odsri] is the same as [mstart]
c
c     work(6): [prt] Pointer of parent of current node. e.g. parent
c           of node specified by pointer "pw" is 
c           work(work(pw + parent))
c           If this value is 0, the current node is root.
c     work(7),work(8): [lc], [rc]: Pointer of left/right child of
c           current node. If both these two value is 0, the current
c           node is leaf.
c           If lc = -1, the left child of current node is
c           not yet specified.
c     work(9), [leaf]: work(pw+leaf) = 0 means not leaf
c           work(pw+leaf) = 1 means leaf
c
c     Next two parameter stores the low rank approximation of
c     off-diagonal block in the same row as node. Here we
c     use ID.
c     work(10): [iL] Pointer of "List"
c     work(11): [iP] Pointer of "Proj"
c           e.g. Proj(work(pw+iP)) is the projection matrix
c           of pnode
c           Proj matrix has the size: p*msize

c
c     work(12) - work(30) is reserved for other attributions.
c
c
c     work(31) - work(50) is other constants.
c
c     work(31): [kappa]: finest(max) level. Given by
c           kappa = INT(log_2(n/2/p))
c           In practice, kappa = INT(log(n*done/2/p)/log(2.0)+eps)
c     work(32): [sns]: Single Node Size. Number of memories need to
c           store a single node.
c           sns = number of memories used from work(2) to work(30)
c     work(33): [WE]: End pointer for "work".
c     work(34): [LE]: End pointer for "List"
c     work(35): [PE]: End pointer for "Proj"
c     work(36): n
c     work(37): p
c
c
      
      integer List(n*n)
c     List: Record A's columns that is selected in ID
      complex *16 Proj(n*n)
c     Proj: Record projection matrices in ID in compressed form
c           (does not contains indentity matrices).
c
c
      call CPU_TIME(ts)
      TIME(1) = ts
      call constructTree(A,LDA,work,List,Proj,n,p)
c     call prinf('List = *',List,work(34)-1)
c     call prin2('Proj = *',Proj,(work(35)-1)*2)
      call CPU_TIME(te)
      TIME(2) = te - ts
      

cc    call prinf('work = *', work(51), work(33)-51)

      call CPU_TIME(ts)
      call Det(A,LDA,work,List,Proj,de)
      call CPU_TIME(te)
      TIME(3) = te - ts
      

      TIME(1) = te - TIME(1)
      


      end
      

      subroutine Det(A,LDA,work,List,Proj,de)
      integer LDA,List(*),work(*)
      complex *16 A(LDA,*),Proj(*)
      complex *16 b(LDA)
      real *8 de

      integer i,j,k,l,pw1,pw2,irow,nrow,pl,pp
      integer kappa, sns, n, p
      integer, allocatable :: IPIV(:)
      complex *16, allocatable :: LU(:,:)

      de = 1.0
      kappa = work(31)
      sns = work(32)
      n = work(36)
      p = work(37)
      pw1 = (2**kappa-1)*sns+work(1) !the first leaf
      allocate(IPIV(2*work(pw1+work(4))))
      allocate(LU(LDA,2*work(pw1+work(4))))


      
      

c     Diagonal Part
      do k = 1, 2**kappa
c           call prinf('kappa = *',kappa,1)
            pw2 = pw1
            irow = work(pw1 + work(3))
            nrow = work(pw1 + work(4))

c     Computing Det of diag submatrix using LU
            call ZGETRF_2(nrow,a(irow,irow),LDA,LU,IPIV)
            
            do j = 1,nrow
                  de = de*LU(j,j)
                  if (IPIV(j) /= j) then
                        de = - de
                  endif
            enddo !j
c           call prin2('det = *',de,1)

c     Update off-diag part of A
            do j = kappa,1,-1
                  pl = work(pw2 + work(10))
c                 call prinf('pl = *',pl,1)
                  do i = 1,p
                        icol = List(pl+i-1)
c                       call prinf('irow = *',irow,1)
c                       call prinf('icol = *',icol,1)
c                       call prinf('nrow = *',nrow,1)
                        call LUBackSolveZ(nrow,1,LU,
     &                  LDA,IPIV,A(irow,icol),LDA)
                  enddo !i

                  pw2 = work(pw2 + work(6))
            enddo !j

            pw1 = pw1 + sns
      enddo !k
c     call prin2('A = *',A,LDA*LDA*2)
      
      deallocate(IPIV)
      deallocate(LU)
      
c     Off-diag low-rank part
      do l = kappa-1,0,-1
            pw1 = (2**l-1)*sns+work(1) !first node in level l
            allocate(IPIV(2*work(pw1 + work(4))))
            allocate(LU(LDA,work(pw1 + work(4))))
            
            do k = 1,2**l
                  pw2 = pw1
                  irow = work(pw1 + work(3))
                  nrow = work(pw1 + work(4))

c     Compute determinant by sylvester's determinant theorem                  
                  call sylvDet_ID(A,LDA,work,List,Proj,de,pw1)

c     Update off-diagonal low rank matrix
                  do j = l,1,-1
                        pl = work(pw2 + work(10))
                        do i = 1,p
                              icol = List(pl+i-1)
                              call SMW_ID(A,LDA,work,List,Proj,
     &                        A(irow,icol),LDA,1,pw1)
                        enddo !i

                        pw2 = work(pw2 + work(6))

                  enddo !j

                  pw1 = pw1 + sns

            enddo !k
            deallocate(IPIV)
            deallocate(LU)
      enddo !l

      


      return
      end

      subroutine sylvDet_ID(A,LDA,work,List,Proj,de,pw)
      implicit none
      integer LDA, List(*), work(*),pw
      complex *16 A(LDA,*), Proj(*)
      real *8 de

      complex *16,allocatable :: U1(:,:),U2(:,:),VT1(:,:),VT2(:,:)
      complex *16,allocatable :: prod1(:,:),prod2(:,:),B(:,:)
      integer INFO, i,j,k
      integer p,icol,irow,nL,nR,pl,pp,pwL,pwR
      integer,allocatable :: IPIV(:)


      p = work(37)
      pwL = work(pw + work(7)) ! lc of pw
      pwR = work(pw + work(8)) ! rc of pw
      nL = work(pwL + work(4))
      nR = work(pwR + work(4))

      allocate(U1(nL,p))
      allocate(U2(nR,p))
      allocate(VT1(p,nR))
      allocate(VT2(p,nL))
      allocate(prod1(p,p))
      allocate(prod2(p,p))
      allocate(B(2*p,2*p))
      allocate(IPIV(2*p))

c     call prinf('p = *',p,1)
c     call prinf('nL = *',nL,1)
c     call prinf('nR = *',nR,1)

c     Construct upper-left block of U (copies of submatrix of A)
      pl = work(pwL+work(10))
c     call prinf('pl = *',pl,1)
      do j = 1,p
            icol = List(pl+j-1)
            do i = 1,nL
                  irow = work(pwL+work(3))+i-1
c                 call prinf('irow = *',irow,1)
c                 call prinf('icol = *',icol,1)
                  U1(i,j)=a(irow,icol)
            enddo !i
      enddo !j
      
c     Construct lower-right block of U
      pl = work(pwR+work(10))
      do j = 1,p
            icol = List(pl+j-1)
            do i = 1,nR
                  irow = work(pwR+work(3))+i-1
                  U2(i,j)=a(irow,icol)
            enddo !i
      enddo !j

c     Construct upper-right block of VT
      pp = work(pwL+work(11))
      do j=1,nR
            do i=1,p
                  VT1(i,j) = Proj(pp+i+(j-1)*p-1)
            enddo
      enddo

c     Construct lower-left block of VT
      pp = work(pwR+work(11))
      do j=1,nL
            do i=1,p
                  VT2(i,j) = Proj(pp+i+(j-1)*p-1)
            enddo
      enddo

      prod1 = matmul(VT1,U2)
      prod2 = matmul(VT2,U1)

c     call prin2('U1 = *',U1,p*nL*2)
c     call prin2('U2 = *',U2,p*nR*2)
c     call prin2('VT1 = *',VT1,p*nR*2)
c     call prin2('VT2 = *',VT2,p*nL*2)
c     call prin2('prod1 = *',prod1,p*p*2)
c     call prin2('prod2 = *',prod2,p*p*2)

      do i = 1,p
            do j = 1,p
                  B(i,j+p) = prod1(i,j)
                  B(i+p,j) = prod2(i,j)
                  if (i == j) then
                        B(i,j) = 1
                  else
                        B(i,j) = 0
                  endif
            enddo !j
      enddo !i

      do i = p+1,p*2
            do j = p+1,p*2
                  if (i == j) then
                        B(i,j) = 1
                  else
                        B(i,j) = 0
                  endif
            enddo !j
      enddo !i

c     call prin2('B = *',B,p*p*8)

      call ZGETRF(p*2,p*2,B,p*2,IPIV,INFO)

      do i = 1,2*p
            de = de * B(i,i)
            if (IPIV(i) /= i) then
                  de = - de
            endif

      enddo

      deallocate(U1)
      deallocate(U2)
      deallocate(VT1)
      deallocate(VT2)
      deallocate(prod1)
      deallocate(prod2)
      deallocate(B)
      deallocate(IPIV)

      return
      end
