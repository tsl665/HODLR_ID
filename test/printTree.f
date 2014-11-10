      subroutine printTree(node, sns, tnn)
c     node: tree structure
c     sns: single node size
c     tnn: total node number
      integer node(1)
      integer sns, tnn
      integer i, j, k

      integer nnum
      
      do k = 1, tnn
            nnum = k - 1
            print *, 'Node', nnum, ': '
c           do j = 0, sns - 1
c                 i = (k-1)*sns+1+j
c                 print *, i, ': ', node(i)
c           enddo
            j = nnum*sns + 1
            print *, 'Level:  ',j, node(j + 0)
            print *, 'mstart: ',j+1 ,node(j + 1)
            print *, 'msize:  ',j+2, node(j + 2)
            print *, 'lr:     ',j+3, node(j + 3)
            print *, 'odsri:  ',j+4, node(j + 4)
            print *, 'odsci:  ',j+5, node(j + 5)
            print *, 'prt:    ',j+6, node(j + 6)
            print *, 'lc:     ',j+7, node(j + 7)
            print *, 'rc:     ',j+8, node(j + 8)
            print *, 'leaf:   ',j+9, node(j + 9)
            print *, 'iP:     ',j+10, node(j + 10)
            print *, 'iL:     ',j+11, node(j + 11)

            print *, '  '
      enddo

      return
      end
