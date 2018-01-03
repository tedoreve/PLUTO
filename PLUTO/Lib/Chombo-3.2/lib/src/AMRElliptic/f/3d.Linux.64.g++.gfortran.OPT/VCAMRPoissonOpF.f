      subroutine GSRBHELMHOLTZVC3D(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,nphicomp
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,nrhscomp
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,dx
     & ,alpha
     & ,aCoef
     & ,iaCoeflo0,iaCoeflo1,iaCoeflo2
     & ,iaCoefhi0,iaCoefhi1,iaCoefhi2
     & ,naCoefcomp
     & ,beta
     & ,bCoef0
     & ,ibCoef0lo0,ibCoef0lo1,ibCoef0lo2
     & ,ibCoef0hi0,ibCoef0hi1,ibCoef0hi2
     & ,nbCoef0comp
     & ,bCoef1
     & ,ibCoef1lo0,ibCoef1lo1,ibCoef1lo2
     & ,ibCoef1hi0,ibCoef1hi1,ibCoef1hi2
     & ,nbCoef1comp
     & ,bCoef2
     & ,ibCoef2lo0,ibCoef2lo1,ibCoef2lo2
     & ,ibCoef2hi0,ibCoef2hi1,ibCoef2hi2
     & ,nbCoef2comp
     & ,lambda
     & ,ilambdalo0,ilambdalo1,ilambdalo2
     & ,ilambdahi0,ilambdahi1,ilambdahi2
     & ,nlambdacomp
     & ,redBlack
     & )
      implicit none
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & 0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & irhslo2:irhshi2,
     & 0:nrhscomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL*8 dx
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0,iaCoeflo1,iaCoeflo2
      integer iaCoefhi0,iaCoefhi1,iaCoefhi2
      REAL*8 aCoef(
     & iaCoeflo0:iaCoefhi0,
     & iaCoeflo1:iaCoefhi1,
     & iaCoeflo2:iaCoefhi2,
     & 0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1,ibCoef0lo2
      integer ibCoef0hi0,ibCoef0hi1,ibCoef0hi2
      REAL*8 bCoef0(
     & ibCoef0lo0:ibCoef0hi0,
     & ibCoef0lo1:ibCoef0hi1,
     & ibCoef0lo2:ibCoef0hi2,
     & 0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1,ibCoef1lo2
      integer ibCoef1hi0,ibCoef1hi1,ibCoef1hi2
      REAL*8 bCoef1(
     & ibCoef1lo0:ibCoef1hi0,
     & ibCoef1lo1:ibCoef1hi1,
     & ibCoef1lo2:ibCoef1hi2,
     & 0:nbCoef1comp-1)
      integer nbCoef2comp
      integer ibCoef2lo0,ibCoef2lo1,ibCoef2lo2
      integer ibCoef2hi0,ibCoef2hi1,ibCoef2hi2
      REAL*8 bCoef2(
     & ibCoef2lo0:ibCoef2hi0,
     & ibCoef2lo1:ibCoef2hi1,
     & ibCoef2lo2:ibCoef2hi2,
     & 0:nbCoef2comp-1)
      integer nlambdacomp
      integer ilambdalo0,ilambdalo1,ilambdalo2
      integer ilambdahi0,ilambdahi1,ilambdahi2
      REAL*8 lambda(
     & ilambdalo0:ilambdahi0,
     & ilambdalo1:ilambdahi1,
     & ilambdalo2:ilambdahi2,
     & 0:nlambdacomp-1)
      integer redBlack
      REAL*8 dxinv,lofphi
      integer n,ncomp,indtot,imin,imax
      integer i,j,k
      ncomp = nphicomp
      if (ncomp .ne. nphicomp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef1comp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef2comp) then
         call MAYDAYERROR()
      endif
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp - 1
        do k=iregionlo2, iregionhi2
          do j=iregionlo1, iregionhi1
            imin = iregionlo0
            indtot = imin + j + k
            imin = imin + abs(mod(indtot + redBlack, 2))
            imax = iregionhi0
            do i = imin, imax, 2
              lofphi =
     & alpha * aCoef(i,j,k,n) * phi(i,j,k,n)
     & - beta *
     & (
     & bCoef0(i+1,j ,k ,n)
     & * (phi(i+1,j ,k ,n)-phi(i ,j ,k ,n))
     & - bCoef0(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n)-phi(i-1,j ,k ,n))
     & + bCoef1(i ,j+1,k ,n)
     & * (phi(i ,j+1,k ,n)-phi(i ,j ,k ,n))
     & - bCoef1(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n)-phi(i ,j-1,k ,n))
     & + bCoef2(i ,j ,k+1,n)
     & * (phi(i ,j ,k+1,n)-phi(i ,j ,k ,n))
     & - bCoef2(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n)-phi(i ,j ,k-1,n))
     & ) * dxinv
              phi(i,j,k,n) = phi(i,j,k,n)
     & - lambda(i,j,k,n) * (lofphi - rhs(i,j,k,n))
            enddo
          enddo
        enddo
      enddo
      return
      end
      subroutine VCCOMPUTEOP3D(
     & lofphi
     & ,ilofphilo0,ilofphilo1,ilofphilo2
     & ,ilofphihi0,ilofphihi1,ilofphihi2
     & ,nlofphicomp
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,nphicomp
     & ,alpha
     & ,aCoef
     & ,iaCoeflo0,iaCoeflo1,iaCoeflo2
     & ,iaCoefhi0,iaCoefhi1,iaCoefhi2
     & ,naCoefcomp
     & ,beta
     & ,bCoef0
     & ,ibCoef0lo0,ibCoef0lo1,ibCoef0lo2
     & ,ibCoef0hi0,ibCoef0hi1,ibCoef0hi2
     & ,nbCoef0comp
     & ,bCoef1
     & ,ibCoef1lo0,ibCoef1lo1,ibCoef1lo2
     & ,ibCoef1hi0,ibCoef1hi1,ibCoef1hi2
     & ,nbCoef1comp
     & ,bCoef2
     & ,ibCoef2lo0,ibCoef2lo1,ibCoef2lo2
     & ,ibCoef2hi0,ibCoef2hi1,ibCoef2hi2
     & ,nbCoef2comp
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,dx
     & )
      implicit none
      integer nlofphicomp
      integer ilofphilo0,ilofphilo1,ilofphilo2
      integer ilofphihi0,ilofphihi1,ilofphihi2
      REAL*8 lofphi(
     & ilofphilo0:ilofphihi0,
     & ilofphilo1:ilofphihi1,
     & ilofphilo2:ilofphihi2,
     & 0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & 0:nphicomp-1)
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0,iaCoeflo1,iaCoeflo2
      integer iaCoefhi0,iaCoefhi1,iaCoefhi2
      REAL*8 aCoef(
     & iaCoeflo0:iaCoefhi0,
     & iaCoeflo1:iaCoefhi1,
     & iaCoeflo2:iaCoefhi2,
     & 0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1,ibCoef0lo2
      integer ibCoef0hi0,ibCoef0hi1,ibCoef0hi2
      REAL*8 bCoef0(
     & ibCoef0lo0:ibCoef0hi0,
     & ibCoef0lo1:ibCoef0hi1,
     & ibCoef0lo2:ibCoef0hi2,
     & 0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1,ibCoef1lo2
      integer ibCoef1hi0,ibCoef1hi1,ibCoef1hi2
      REAL*8 bCoef1(
     & ibCoef1lo0:ibCoef1hi0,
     & ibCoef1lo1:ibCoef1hi1,
     & ibCoef1lo2:ibCoef1hi2,
     & 0:nbCoef1comp-1)
      integer nbCoef2comp
      integer ibCoef2lo0,ibCoef2lo1,ibCoef2lo2
      integer ibCoef2hi0,ibCoef2hi1,ibCoef2hi2
      REAL*8 bCoef2(
     & ibCoef2lo0:ibCoef2hi0,
     & ibCoef2lo1:ibCoef2hi1,
     & ibCoef2lo2:ibCoef2hi2,
     & 0:nbCoef2comp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL*8 dx
      REAL*8 dxinv
      integer n,ncomp
      integer i,j,k
      ncomp = nphicomp
      if (ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef1comp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef2comp) then
         call MAYDAYERROR()
      endif
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          lofphi(i,j,k,n) =
     & alpha * aCoef(i,j,k,n) * phi(i,j,k,n)
     & - beta *
     & (
     & bCoef0(i+1,j ,k ,n)
     & * (phi(i+1,j ,k ,n) - phi(i ,j ,k ,n))
     & - bCoef0(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n) - phi(i-1,j ,k ,n))
     & + bCoef1(i ,j+1,k ,n)
     & * (phi(i ,j+1,k ,n) - phi(i ,j ,k ,n))
     & - bCoef1(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n) - phi(i ,j-1,k ,n))
     & + bCoef2(i ,j ,k+1,n)
     & * (phi(i ,j ,k+1,n) - phi(i ,j ,k ,n))
     & - bCoef2(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n) - phi(i ,j ,k-1,n))
     & ) * dxinv
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine VCCOMPUTERES3D(
     & res
     & ,ireslo0,ireslo1,ireslo2
     & ,ireshi0,ireshi1,ireshi2
     & ,nrescomp
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,nphicomp
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,nrhscomp
     & ,alpha
     & ,aCoef
     & ,iaCoeflo0,iaCoeflo1,iaCoeflo2
     & ,iaCoefhi0,iaCoefhi1,iaCoefhi2
     & ,naCoefcomp
     & ,beta
     & ,bCoef0
     & ,ibCoef0lo0,ibCoef0lo1,ibCoef0lo2
     & ,ibCoef0hi0,ibCoef0hi1,ibCoef0hi2
     & ,nbCoef0comp
     & ,bCoef1
     & ,ibCoef1lo0,ibCoef1lo1,ibCoef1lo2
     & ,ibCoef1hi0,ibCoef1hi1,ibCoef1hi2
     & ,nbCoef1comp
     & ,bCoef2
     & ,ibCoef2lo0,ibCoef2lo1,ibCoef2lo2
     & ,ibCoef2hi0,ibCoef2hi1,ibCoef2hi2
     & ,nbCoef2comp
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,dx
     & )
      implicit none
      integer nrescomp
      integer ireslo0,ireslo1,ireslo2
      integer ireshi0,ireshi1,ireshi2
      REAL*8 res(
     & ireslo0:ireshi0,
     & ireslo1:ireshi1,
     & ireslo2:ireshi2,
     & 0:nrescomp-1)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & 0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & irhslo2:irhshi2,
     & 0:nrhscomp-1)
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0,iaCoeflo1,iaCoeflo2
      integer iaCoefhi0,iaCoefhi1,iaCoefhi2
      REAL*8 aCoef(
     & iaCoeflo0:iaCoefhi0,
     & iaCoeflo1:iaCoefhi1,
     & iaCoeflo2:iaCoefhi2,
     & 0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1,ibCoef0lo2
      integer ibCoef0hi0,ibCoef0hi1,ibCoef0hi2
      REAL*8 bCoef0(
     & ibCoef0lo0:ibCoef0hi0,
     & ibCoef0lo1:ibCoef0hi1,
     & ibCoef0lo2:ibCoef0hi2,
     & 0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1,ibCoef1lo2
      integer ibCoef1hi0,ibCoef1hi1,ibCoef1hi2
      REAL*8 bCoef1(
     & ibCoef1lo0:ibCoef1hi0,
     & ibCoef1lo1:ibCoef1hi1,
     & ibCoef1lo2:ibCoef1hi2,
     & 0:nbCoef1comp-1)
      integer nbCoef2comp
      integer ibCoef2lo0,ibCoef2lo1,ibCoef2lo2
      integer ibCoef2hi0,ibCoef2hi1,ibCoef2hi2
      REAL*8 bCoef2(
     & ibCoef2lo0:ibCoef2hi0,
     & ibCoef2lo1:ibCoef2hi1,
     & ibCoef2lo2:ibCoef2hi2,
     & 0:nbCoef2comp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL*8 dx
      REAL*8 dxinv
      integer n,ncomp
      integer i,j,k
      ncomp = nphicomp
      if (ncomp .ne. nrescomp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef0comp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef1comp) then
         call MAYDAYERROR()
      endif
      if (ncomp .ne. nbCoef2comp) then
         call MAYDAYERROR()
      endif
      dxinv = (1.0d0)/(dx*dx)
      do n = 0, ncomp-1
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          res(i,j,k,n) =
     & rhs(i,j,k,n)
     & - (alpha * aCoef(i,j,k,n) * phi(i,j,k,n)
     & - beta *
     & (
     & bCoef0(i+1,j ,k ,n)
     & * (phi(i+1,j ,k ,n) - phi(i ,j ,k ,n))
     & - bCoef0(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n) - phi(i-1,j ,k ,n))
     & + bCoef1(i ,j+1,k ,n)
     & * (phi(i ,j+1,k ,n) - phi(i ,j ,k ,n))
     & - bCoef1(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n) - phi(i ,j-1,k ,n))
     & + bCoef2(i ,j ,k+1,n)
     & * (phi(i ,j ,k+1,n) - phi(i ,j ,k ,n))
     & - bCoef2(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n) - phi(i ,j ,k-1,n))
     & ) * dxinv
     & )
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine RESTRICTRESVC3D(
     & res
     & ,ireslo0,ireslo1,ireslo2
     & ,ireshi0,ireshi1,ireshi2
     & ,nrescomp
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,nphicomp
     & ,rhs
     & ,irhslo0,irhslo1,irhslo2
     & ,irhshi0,irhshi1,irhshi2
     & ,nrhscomp
     & ,alpha
     & ,aCoef
     & ,iaCoeflo0,iaCoeflo1,iaCoeflo2
     & ,iaCoefhi0,iaCoefhi1,iaCoefhi2
     & ,naCoefcomp
     & ,beta
     & ,bCoef0
     & ,ibCoef0lo0,ibCoef0lo1,ibCoef0lo2
     & ,ibCoef0hi0,ibCoef0hi1,ibCoef0hi2
     & ,nbCoef0comp
     & ,bCoef1
     & ,ibCoef1lo0,ibCoef1lo1,ibCoef1lo2
     & ,ibCoef1hi0,ibCoef1hi1,ibCoef1hi2
     & ,nbCoef1comp
     & ,bCoef2
     & ,ibCoef2lo0,ibCoef2lo1,ibCoef2lo2
     & ,ibCoef2hi0,ibCoef2hi1,ibCoef2hi2
     & ,nbCoef2comp
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,dx
     & )
      implicit none
      integer nrescomp
      integer ireslo0,ireslo1,ireslo2
      integer ireshi0,ireshi1,ireshi2
      REAL*8 res(
     & ireslo0:ireshi0,
     & ireslo1:ireshi1,
     & ireslo2:ireshi2,
     & 0:nrescomp-1)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & 0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL*8 rhs(
     & irhslo0:irhshi0,
     & irhslo1:irhshi1,
     & irhslo2:irhshi2,
     & 0:nrhscomp-1)
      REAL*8 alpha
      integer naCoefcomp
      integer iaCoeflo0,iaCoeflo1,iaCoeflo2
      integer iaCoefhi0,iaCoefhi1,iaCoefhi2
      REAL*8 aCoef(
     & iaCoeflo0:iaCoefhi0,
     & iaCoeflo1:iaCoefhi1,
     & iaCoeflo2:iaCoefhi2,
     & 0:naCoefcomp-1)
      REAL*8 beta
      integer nbCoef0comp
      integer ibCoef0lo0,ibCoef0lo1,ibCoef0lo2
      integer ibCoef0hi0,ibCoef0hi1,ibCoef0hi2
      REAL*8 bCoef0(
     & ibCoef0lo0:ibCoef0hi0,
     & ibCoef0lo1:ibCoef0hi1,
     & ibCoef0lo2:ibCoef0hi2,
     & 0:nbCoef0comp-1)
      integer nbCoef1comp
      integer ibCoef1lo0,ibCoef1lo1,ibCoef1lo2
      integer ibCoef1hi0,ibCoef1hi1,ibCoef1hi2
      REAL*8 bCoef1(
     & ibCoef1lo0:ibCoef1hi0,
     & ibCoef1lo1:ibCoef1hi1,
     & ibCoef1lo2:ibCoef1hi2,
     & 0:nbCoef1comp-1)
      integer nbCoef2comp
      integer ibCoef2lo0,ibCoef2lo1,ibCoef2lo2
      integer ibCoef2hi0,ibCoef2hi1,ibCoef2hi2
      REAL*8 bCoef2(
     & ibCoef2lo0:ibCoef2hi0,
     & ibCoef2lo1:ibCoef2hi1,
     & ibCoef2lo2:ibCoef2hi2,
     & 0:nbCoef2comp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL*8 dx
      REAL*8 denom,dxinv,lofphi
      integer n,ncomp
      integer i,j,k
      integer ii,jj,kk
      ncomp = nphicomp
      dxinv = (1.0d0) / (dx*dx)
      denom = 2 *2 *2
      do n = 0, ncomp-1
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0
          ii = i/2
          jj = j/2
          kk = k/2
          lofphi =
     & alpha * aCoef(i,j,k,n) * phi(i,j,k,n)
     & - beta *
     & (
     & bCoef0(i+1,j ,k ,n)
     & * (phi(i+1,j ,k ,n)-phi(i ,j ,k ,n))
     & - bCoef0(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n)-phi(i-1,j ,k ,n))
     & + bCoef1(i ,j+1,k ,n)
     & * (phi(i ,j+1,k ,n)-phi(i ,j ,k ,n))
     & - bCoef1(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n)-phi(i ,j-1,k ,n))
     & + bCoef2(i ,j ,k+1,n)
     & * (phi(i ,j ,k+1,n)-phi(i ,j ,k ,n))
     & - bCoef2(i ,j ,k ,n)
     & * (phi(i ,j ,k ,n)-phi(i ,j ,k-1,n))
     & ) * dxinv
          res(ii,jj,kk,n) = res(ii,jj,kk,n)
     & + (rhs(i,j,k,n) - lofphi) / denom
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine SUMFACES(
     & lhs
     & ,ilhslo0,ilhslo1,ilhslo2
     & ,ilhshi0,ilhshi1,ilhshi2
     & ,nlhscomp
     & ,beta
     & ,bCoefs
     & ,ibCoefslo0,ibCoefslo1,ibCoefslo2
     & ,ibCoefshi0,ibCoefshi1,ibCoefshi2
     & ,nbCoefscomp
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & ,dir
     & ,scale
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nlhscomp
      integer ilhslo0,ilhslo1,ilhslo2
      integer ilhshi0,ilhshi1,ilhshi2
      REAL*8 lhs(
     & ilhslo0:ilhshi0,
     & ilhslo1:ilhshi1,
     & ilhslo2:ilhshi2,
     & 0:nlhscomp-1)
      REAL*8 beta
      integer nbCoefscomp
      integer ibCoefslo0,ibCoefslo1,ibCoefslo2
      integer ibCoefshi0,ibCoefshi1,ibCoefshi2
      REAL*8 bCoefs(
     & ibCoefslo0:ibCoefshi0,
     & ibCoefslo1:ibCoefshi1,
     & ibCoefslo2:ibCoefshi2,
     & 0:nbCoefscomp-1)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer dir
      REAL*8 scale
      REAL*8 sumVal
      integer i,j,k
      integer ii,jj,kk
      integer n
      ii = CHF_ID(0,dir)
      jj = CHF_ID(1,dir)
      kk = CHF_ID(2,dir)
      do n = 0, nlhscomp-1
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0
          sumVal = bCoefs(i+ii,j+jj,k+kk,n)
     & + bCoefs(i ,j ,k ,n)
          lhs(i,j,k,n) = lhs(i,j,k,n) + scale * beta * sumVal
      enddo
      enddo
      enddo
      enddo
      return
      end
