c	V302: 
c	1. Using standard back solve to avoid unnecessary ZGESV
c	2. Add new return variable TIME
c	3. Put small subroutines into a seperate file (SMW_SGESV2.f)

	subroutine HODLR_ID(n,A,LDA,p,b,LDB,r,TIME)
c	n: size of matrix A
c	A: HODLR matrix to solve
c	LDA: Leading Dimension of A
c	p: specific rank
c	b: r.h.s
c	LDB: Leading Dimension of b
c	r: number of r.h.s
c	TIME: an real array that records Elapsed CPU time for each
c		step.
c	TIME(1): total time
c	TIME(2): construct tree
c	TIME(3): dense invert
c	TIME(4): low-rank invert
	implicit none
	integer WLEN
	parameter (WLEN = 1E7) ! Length of work
	integer n, LDA, LDB, p, r
	real *8 TIME(10),ts,te
	complex *16 A(LDA,n)
	complex *16 b(LDB,r)
	
	integer work(WLEN)
c	work: the object stores the tree structure of HODLR matrix.
c	work(1): starting index of tree structure. The first 1 to 
c		work(1)-1 memories are for constants. Here work(1) = 51.
c	
c	work(2) - work(20) is constants about attribution of each
c	single node. Let "tag" be the information we want. "pw"
c	is pointer for current node. Then we access that information
c	by work(pw+tag). "tag" is one of the following:
c
c	work(2): [lvl] Level of current node
c	work(3): [mstart] The starting column(row) of current node
c		in A
c	work(4): [msize] Size of current node
c	work(5): [odsci] Off-Diagonal Starting Column Index
c		Notice that the [odsri] is the same as [mstart]
c
c	work(6): [prt] Pointer of parent of current node. e.g. parent
c		of node specified by pointer "pw" is 
c		work(work(pw + parent))
c		If this value is 0, the current node is root.
c	work(7),work(8): [lc], [rc]: Pointer of left/right child of
c		current node. If both these two value is 0, the current
c		node is leaf.
c		If lc = -1, the left child of current node is
c		not yet specified.
c	work(9), [leaf]: work(pw+leaf) = 0 means not leaf
c		work(pw+leaf) = 1 means leaf
c
c	Next two parameter stores the low rank approximation of
c	off-diagonal block in the same row as node. Here we
c	use ID.
c	work(10): [iL] Pointer of "List"
c	work(11): [iP] Pointer of "Proj"
c		e.g. Proj(work(pw+iP)) is the projection matrix
c		of pnode
c		Proj matrix has the size: p*msize

c
c	work(12) - work(30) is reserved for other attributions.
c
c
c	work(31) - work(50) is other constants.
c
c	work(31): [kappa]: finest(max) level. Given by
c		kappa = INT(log_2(n/2/p))
c		In practice, kappa = INT(log(n*done/2/p)/log(2.0)+eps)
c	work(32): [sns]: Single Node Size. Number of memories need to
c		store a single node.
c		sns = number of memories used from work(2) to work(30)
c	work(33): [WE]: End pointer for "work".
c	work(34): [LE]: End pointer for "List"
c	work(35): [PE]: End pointer for "Proj"
c	work(36): n
c	work(37): p
c
c
	
	integer List(n*n)
c	List: Record A's columns that is selected in ID
	complex *16 Proj(n*n)
c	Proj: Record projection matrices in ID in compressed form
c		(does not contains indentity matrices).
c
c
	call CPU_TIME(ts)
	TIME(1) = ts
	call constructTree(A,LDA,work,List,Proj,n,p)
	call CPU_TIME(te)
	TIME(2) = te - ts
	

cc	call prinf('work = *', work(51), work(33)-51)

	call CPU_TIME(ts)
	call DenseInverse(A,LDA,work,List,Proj,B,LDB,r)
	call CPU_TIME(te)
	TIME(3) = te - ts
	

	call CPU_TIME(ts)
	call lrInverse(A,LDA,work,List,Proj,B,LDB,r)
	call CPU_TIME(te)
	TIME(4) = te - ts
	TIME(1) = te - TIME(1)
	


	end

	subroutine constructTree(A,LDA,work,List,Proj,n,p)
	implicit none
	integer work(*),List(*),LDA,n,p
	complex *16 Proj(*), A(LDA,*)

	real *8 done,eps
	complex *16 ima
	parameter (done = 1.0, ima = (0,1) ,eps = 1E-7)	


	integer i,j,k,root
	integer pp,pl,pw

c	Initialization
	work(1) = 51
	do i = 2, 11
		work(i) = i-2
	enddo
	work(31) = INT (log(n*done/2/p)/log(2.0)+eps) !kappa
	work(32) = 10          !sns
	work(33) = work(1) - 1 !WE
	work(34) = 1           !LE
	work(35) = 1           !PE
	work(36) = n
	work(37) = p

c	Construct root node
	root = work(1)
	work(root + work(2)) = 0  !lvl
	work(root + work(3)) = 1  !mstart
	work(root + work(4)) = n  !msize
	work(root + work(5)) = 0  !odsci(no low-rank part)
	work(root + work(6)) = 0  !prt(no parent)
	work(root + work(7)) = -1 !lc
	work(root + work(8)) = -1 !rc
	work(root + work(9)) = 0  !leaf
	work(root + work(10)) = 0 !iL
	work(root + work(11)) = 0 !iP

	work(33) = root + work(32) ! WE = root + sns

	pw = root
cc	pp = 1
cc	pl = 1

c	Construct the other nodes
	do
		if (pw >= work(33)) exit

		if (work(pw+work(7)) == -1) then
			call newNode(A,LDA,work, List, Proj, pw, 'L')
		endif

		if (work(pw+work(8)) == -1) then
			call newNode(A, LDA, work, List, Proj, pw, 'R')
		endif

		pw = pw + work(32) ! pw = pw + sns

	enddo

	return
	end

	subroutine newNode(A, LDA, work,List,Proj,pw,child)
	implicit none
	integer LDA, pw, work(*),List(*)
	complex *16 A(LDA, *), Proj(*)
	character child
	
	integer pw1, pw2, pl, pp, i,j,k, nL, nR, nS, irow, icol
	integer kappa, n, p

	pw1 = pw
	pw2 = work(33) !WE
	pl = work(34)
	pp = work(35)
	kappa = work(31)
	n = work(36)
	p = work(37)
	nL = INT(work(pw1 + work(4))/2) ! left child size
	nR = work(pw1 + work(4)) - nL   ! right child size
	nS = work(pw1 + work(3))        ! parent starting index

c	Update independent attribution
	work(pw2 + work(2)) = work(pw1 + work(2)) + 1 !lvl
	work(pw2 + work(6)) = pw1 !prt
	work(pw2 + work(10)) = pl
	work(pw2 + work(11)) = pp

c	Update attribution depending on level
	if (work(pw2 + work(2)) < kappa) then
		work(pw2 + work(7)) = -1 !lc
		work(pw2 + work(8)) = -1 !rc
		work(pw2 + work(9)) = 0  !leaf
	else
		work(pw2 + work(7)) = 0
		work(pw2 + work(8)) = 0
		work(pw2 + work(9)) = 1
	endif
	
c	Update attribution depending on child type
	if ( child == 'L' ) then
		work(pw2 + work(3)) = nS      !mstart
		work(pw2 + work(4)) = nL      !msize
		work(pw2 + work(5)) = nS + nL !odsci
	
		work(pw1 + work(7)) = pw2	!lc of parent
	elseif (child == 'R') then
		work(pw2 + work(3)) = nS + nL !mstart
		work(pw2 + work(4)) = nR	!msize
		work(pw2 + work(5)) = nS	!odsci

		work(pw1 + work(8)) = pw2	!rc of parent
	else
		print *, 'Wrong input "child": Must be "L" or "R"!'
	endif

c	compute low-rank approximation
	irow = work(pw2 + work(3)) ! odsri = mstart
	icol = work(pw2 + work(5)) ! odsci
	if ( child == 'L') then
		call IDZR_ID_2(nL,nR,A(irow,icol),LDA,p,icol,
     &	List(pl),Proj(pp))
	elseif (child == 'R') then
		call IDZR_ID_2(nR,nL,A(irow,icol),LDA,p,icol,
     &	List(pl),Proj(pp))
	endif


c	Update pointers
	work(34) = work(34) + p !LE

	if (child == 'L') then  !PE
		work(35) = work(35) + nR*p
	elseif (child == 'R') then
		work(35) = work(35) + nL*p
	endif

	work(33) = work(33) + work(32) !WE = WE + sns

	return

	end

	subroutine IDZR_ID_2(m,n,A,LDA,krank,lrcol,List,Proj)
	implicit none
	integer LDA, List(*), Ltemp(n),krank,lrcol,m,n
	complex *16 A(LDA, *), Proj(krank,n)

	real *8 rnorms(n)
	complex *16 Atemp(m,n)
	integer i,j,k

	do j = 1,n
		do i = 1,m
			Atemp(i,j) = A(i,j)
		enddo
	enddo

	call idzr_id(m,n,Atemp,krank,Ltemp,rnorms)

	call recordProj(m,n,krank,Proj,Atemp,Ltemp)

c	record list and adjust to index of A
	do j = 1,krank
		List(j) = Ltemp(j) + lrcol -1
	enddo

c	call checkLR(m,n,krank,A,LDA,Ltemp,Proj)


	return
	end

	subroutine recordProj(m,n,krank,Proj,Atemp,Ltemp)
	integer m,n,krank,Ltemp(*)
	complex *16 Atemp(krank,n-krank), Proj(krank,n)
	integer i,j,icol

c	Identity part:
	do j = 1,krank
		icol = Ltemp(j)
		do i = 1,krank
			if (i == j) then
				proj(i,icol) = 1
			else
				proj(i,icol) = 0
			endif
		enddo !i
	enddo !j

c	The rest of Proj
	do j = krank + 1,n
		icol = Ltemp(j)
		do i = 1,krank
			proj(i,icol) = Atemp(i,j-krank)
		enddo !i
	enddo !j

	return
	end



	subroutine denseInverse(A,LDA,work,List,Proj,B,LDB,r)
	integer LDA,LDB,r, List(*),work(*)
	complex *16 A(LDA,*),Proj(*),B(LDB,r)

	integer i,j,k,pw1,pw2,irow,nrow,pl,pp
	integer kappa, sns, n, p
	integer, allocatable :: IPIV(:)
	complex *16, allocatable :: LU(:,:)

c	integer IPIV(LDA)
c	complex *16 LU(LDA,LDA)


	kappa = work(31)
	sns = work(32)
	n = work(36)
	p = work(37)
	
	pw1 = (2**kappa-1)*sns+work(1) !the first leaf
	allocate(IPIV(2*work(pw1+work(4))))
	allocate(LU(LDA,2*work(pw1+work(4))))

	do k = 1, 2**kappa
		pw2 = pw1
		irow = work(pw1 + work(3))
		nrow = work(pw1 + work(4))
c		allocate(IPIV(2*nrow))
c		allocate(LU(2*nrow,2*nrow))
cccc	check for b: n-by-r
		call ZGESV_3(nrow,r,a(irow,irow),LDA,b(irow,1),
     &	LDB,LU,IPIV)
		
		do j = kappa,1,-1
			pl = work(pw2 + work(10))
			do i = 1,p
				icol = List(pl+i-1)
c				call ZGESV_2(nrow,1,A(irow,irow),LDA,
c     &			a(irow,icol),LDA)
				call LUBackSolveZ(nrow,1,LU,
     &			LDA,IPIV,A(irow,icol),LDA)
			enddo !i

			pw2 = work(pw2 + work(6))
		enddo !j

		pw1 = pw1 + sns
c		deallocate(IPIV)
c		deallocate(LU)
	enddo !k

	deallocate(IPIV)
	deallocate(LU)

	return
	end

	subroutine lrInverse(A,LDA,work,List,Proj,B,LDB,r)
	implicit none
	integer LDA,LDB,r, List(*),work(*)
	complex *16 A(LDA,*),Proj(*),B(LDB,r)

	integer i,j,k,l,pw1,pw2,irow,nrow,icol,pl,pp
	integer kappa, sns, n, p

	kappa = work(31)
	sns = work(32)
	n = work(36)
	p = work(37)

c	call prinf('kappa = *',kappa,1)

	do l = kappa-1,0,-1
		pw1 = (2**l-1)*sns+work(1) !first node in level l
		
		do k = 1,2**l
			pw2 = pw1
			irow = work(pw1 + work(3))
			nrow = work(pw1 + work(4))

			do j = l,1,-1
				pl = work(pw2 + work(10))
				do i = 1,p
					icol = List(pl+i-1)
					call SMW_ID(A,LDA,work,List,Proj,
     &				A(irow,icol),LDA,1,pw1)
				enddo !i

				pw2 = work(pw2 + work(6))

			enddo !j
			call SMW_ID(A,LDA,work,List,Proj,B(irow,1),LDB,r,
     &		pw1)

			pw1 = pw1 + sns

		enddo !k
	enddo !l

	return
	end
