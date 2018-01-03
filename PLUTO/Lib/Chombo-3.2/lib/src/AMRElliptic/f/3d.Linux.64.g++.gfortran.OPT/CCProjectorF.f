      subroutine CCPEXTRAPTODOMFACE(
     & facevel
     & ,ifacevello0,ifacevello1,ifacevello2
     & ,ifacevelhi0,ifacevelhi1,ifacevelhi2
     & ,facedir
     & ,ishift
     & ,ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
     & ,ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifacevello0,ifacevello1,ifacevello2
      integer ifacevelhi0,ifacevelhi1,ifacevelhi2
      REAL*8 facevel(
     & ifacevello0:ifacevelhi0,
     & ifacevello1:ifacevelhi1,
     & ifacevello2:ifacevelhi2)
      integer facedir
      integer ishift
      integer ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
      integer ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
      integer i,j,k
      integer ioff,joff,koff
      ioff = ishift*chf_id(0,facedir)
      joff = ishift*chf_id(1,facedir)
      koff = ishift*chf_id(2,facedir)
      do k = ifaceboxlo2,ifaceboxhi2
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
      facevel(i,j,k) =
     & ((2.0d0)*facevel(i+ ioff,j+ joff,k+ koff)
     & - facevel(i+2*ioff,j+2*joff,k+2*koff))
      enddo
      enddo
      enddo
      return
      end
      subroutine CCPAVECELLTOFACE(
     & facevel
     & ,ifacevello0,ifacevello1,ifacevello2
     & ,ifacevelhi0,ifacevelhi1,ifacevelhi2
     & ,cellvel
     & ,icellvello0,icellvello1,icellvello2
     & ,icellvelhi0,icellvelhi1,icellvelhi2
     & ,facedir
     & ,ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
     & ,ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ifacevello0,ifacevello1,ifacevello2
      integer ifacevelhi0,ifacevelhi1,ifacevelhi2
      REAL*8 facevel(
     & ifacevello0:ifacevelhi0,
     & ifacevello1:ifacevelhi1,
     & ifacevello2:ifacevelhi2)
      integer icellvello0,icellvello1,icellvello2
      integer icellvelhi0,icellvelhi1,icellvelhi2
      REAL*8 cellvel(
     & icellvello0:icellvelhi0,
     & icellvello1:icellvelhi1,
     & icellvello2:icellvelhi2)
      integer facedir
      integer ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
      integer ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
      integer i,j,k
      integer ioff,joff,koff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      do k = ifaceboxlo2,ifaceboxhi2
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
      facevel(i,j,k) =
     & ( cellvel(i ,j ,k )
     & + cellvel(i-ioff,j-joff,k-koff)
     & )*(0.500d0)
      enddo
      enddo
      enddo
      return
      end
      subroutine CCPAVEFACETOCELL(
     & cellgrad
     & ,icellgradlo0,icellgradlo1,icellgradlo2
     & ,icellgradhi0,icellgradhi1,icellgradhi2
     & ,facegrad
     & ,ifacegradlo0,ifacegradlo1,ifacegradlo2
     & ,ifacegradhi0,ifacegradhi1,ifacegradhi2
     & ,facedir
     & ,icellboxlo0,icellboxlo1,icellboxlo2
     & ,icellboxhi0,icellboxhi1,icellboxhi2
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer icellgradlo0,icellgradlo1,icellgradlo2
      integer icellgradhi0,icellgradhi1,icellgradhi2
      REAL*8 cellgrad(
     & icellgradlo0:icellgradhi0,
     & icellgradlo1:icellgradhi1,
     & icellgradlo2:icellgradhi2)
      integer ifacegradlo0,ifacegradlo1,ifacegradlo2
      integer ifacegradhi0,ifacegradhi1,ifacegradhi2
      REAL*8 facegrad(
     & ifacegradlo0:ifacegradhi0,
     & ifacegradlo1:ifacegradhi1,
     & ifacegradlo2:ifacegradhi2)
      integer facedir
      integer icellboxlo0,icellboxlo1,icellboxlo2
      integer icellboxhi0,icellboxhi1,icellboxhi2
      integer i, j, k
      integer ioff, joff, koff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      do k = icellboxlo2,icellboxhi2
      do j = icellboxlo1,icellboxhi1
      do i = icellboxlo0,icellboxhi0
      cellgrad(i,j,k) =
     & (facegrad(i+ioff,j+joff,k+koff)
     & +facegrad(i ,j ,k ))*(0.500d0)
      enddo
      enddo
      enddo
      return
      end
      subroutine CCPCELLGRADFROMFACEDATA(
     & cellgrad
     & ,icellgradlo0,icellgradlo1,icellgradlo2
     & ,icellgradhi0,icellgradhi1,icellgradhi2
     & ,facedata
     & ,ifacedatalo0,ifacedatalo1,ifacedatalo2
     & ,ifacedatahi0,ifacedatahi1,ifacedatahi2
     & ,facedir
     & ,dx
     & ,icellboxlo0,icellboxlo1,icellboxlo2
     & ,icellboxhi0,icellboxhi1,icellboxhi2
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer icellgradlo0,icellgradlo1,icellgradlo2
      integer icellgradhi0,icellgradhi1,icellgradhi2
      REAL*8 cellgrad(
     & icellgradlo0:icellgradhi0,
     & icellgradlo1:icellgradhi1,
     & icellgradlo2:icellgradhi2)
      integer ifacedatalo0,ifacedatalo1,ifacedatalo2
      integer ifacedatahi0,ifacedatahi1,ifacedatahi2
      REAL*8 facedata(
     & ifacedatalo0:ifacedatahi0,
     & ifacedatalo1:ifacedatahi1,
     & ifacedatalo2:ifacedatahi2)
      integer facedir
      REAL*8 dx
      integer icellboxlo0,icellboxlo1,icellboxlo2
      integer icellboxhi0,icellboxhi1,icellboxhi2
      integer i, j, k
      integer ioff, joff, koff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      do k = icellboxlo2,icellboxhi2
      do j = icellboxlo1,icellboxhi1
      do i = icellboxlo0,icellboxhi0
      cellgrad(i,j,k) =
     & (facedata(i+ioff,j+joff,k+koff)
     & -facedata(i ,j ,k ))/dx
      enddo
      enddo
      enddo
      return
      end
