      subroutine HALVEINT(
     & arr
     & ,iarrlo0,iarrlo1,iarrlo2
     & ,iarrhi0,iarrhi1,iarrhi2
     & ,narrcomp
     & ,ibxlo0,ibxlo1,ibxlo2
     & ,ibxhi0,ibxhi1,ibxhi2
     & )
      implicit none
      integer narrcomp
      integer iarrlo0,iarrlo1,iarrlo2
      integer iarrhi0,iarrhi1,iarrhi2
      integer arr(
     & iarrlo0:iarrhi0,
     & iarrlo1:iarrhi1,
     & iarrlo2:iarrhi2,
     & 0:narrcomp-1)
      integer ibxlo0,ibxlo1,ibxlo2
      integer ibxhi0,ibxhi1,ibxhi2
      integer i0,i1,i2
      integer var
      do var = 0, narrcomp-1
      do i2 = ibxlo2,ibxhi2
      do i1 = ibxlo1,ibxhi1
      do i0 = ibxlo0,ibxhi0
            arr(i0,i1,i2, var) = arr(i0,i1,i2, var) / 2
      enddo
      enddo
      enddo
      enddo
      return
      end
