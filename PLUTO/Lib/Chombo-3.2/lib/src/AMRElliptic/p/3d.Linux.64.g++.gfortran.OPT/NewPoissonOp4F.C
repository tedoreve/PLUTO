#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine OPERATORLAP4(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1,ilofphilo2
     &           ,ilofphihi0,ilofphihi1,ilofphihi2
     &           ,nlofphicomp
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,dx
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      integer nlofphicomp
      integer ilofphilo0,ilofphilo1,ilofphilo2
      integer ilofphihi0,ilofphihi1,ilofphihi2
      REAL_T lofphi(
     &           ilofphilo0:ilofphihi0,
     &           ilofphilo1:ilofphihi1,
     &           ilofphilo2:ilofphihi2,
     &           0:nlofphicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      REAL_T dxinv, lap
      integer n,ncomp
      integer i,j,k
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
        
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

          lap = ( 
     &   sixteen*phi(i-1,j  ,k  ,n) - phi(i-2,j  ,k  ,n)
     & + sixteen*phi(i+1,j  ,k  ,n) - phi(i+2,j  ,k  ,n)
     & + sixteen*phi(i  ,j-1,k  ,n) - phi(i  ,j-2,k  ,n)
     & + sixteen*phi(i  ,j+1,k  ,n) - phi(i  ,j+2,k  ,n)
     & + sixteen*phi(i  ,j  ,k-1,n) - phi(i  ,j  ,k-2,n)
     & + sixteen*phi(i  ,j  ,k+1,n) - phi(i  ,j  ,k+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,j,k,n) )
     &       * twelfth * dxinv
          lofphi(i,j,k,n) = alpha*phi(i,j,k,n)+beta*lap
        
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine OPERATORLAPRES4(
     &           r
     &           ,irlo0,irlo1,irlo2
     &           ,irhi0,irhi1,irhi2
     &           ,nrcomp
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2
     &           ,irhshi0,irhshi1,irhshi2
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,dx
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      integer nrcomp
      integer irlo0,irlo1,irlo2
      integer irhi0,irhi1,irhi2
      REAL_T r(
     &           irlo0:irhi0,
     &           irlo1:irhi1,
     &           irlo2:irhi2,
     &           0:nrcomp-1)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      REAL_T dxinv, lap, lhs
      integer n,ncomp
      integer i,j,k
      ncomp = nphicomp
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

          lap = ( 
     &   sixteen*phi(i-1,j  ,k  ,n) - phi(i-2,j  ,k  ,n)
     & + sixteen*phi(i+1,j  ,k  ,n) - phi(i+2,j  ,k  ,n)
     & + sixteen*phi(i  ,j-1,k  ,n) - phi(i  ,j-2,k  ,n)
     & + sixteen*phi(i  ,j+1,k  ,n) - phi(i  ,j+2,k  ,n)
     & + sixteen*phi(i  ,j  ,k-1,n) - phi(i  ,j  ,k-2,n)
     & + sixteen*phi(i  ,j  ,k+1,n) - phi(i  ,j  ,k+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,j,k,n) )
     &       * twelfth * dxinv
          lhs = alpha*phi(i,j,k,n) + beta*lap
          r(i,j,k,n) = rhs(i,j,k,n) - lhs
         
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine RESTRICTRES4(
     &           res
     &           ,ireslo0,ireslo1,ireslo2
     &           ,ireshi0,ireshi1,ireshi2
     &           ,nrescomp
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2
     &           ,irhshi0,irhshi1,irhshi2
     &           ,nrhscomp
     &           ,alpha
     &           ,beta
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,dx
     &           )

      implicit none
      integer nrescomp
      integer ireslo0,ireslo1,ireslo2
      integer ireshi0,ireshi1,ireshi2
      REAL_T res(
     &           ireslo0:ireshi0,
     &           ireslo1:ireshi1,
     &           ireslo2:ireshi2,
     &           0:nrescomp-1)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2,
     &           0:nrhscomp-1)
      REAL_T alpha
      REAL_T beta
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T dx
      REAL_T denom,dxinv,lofphi,lap
      integer n,ncomp
      integer i,j,k
      integer ii,jj,kk
      ncomp = nphicomp
      dxinv = one / (dx*dx)
      denom = D_TERM(2, *2, *2)
      do n = 0, ncomp-1
        
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

          
          ii = i/2 
          jj = j/2 
          kk = k/2 
          lap = ( 
     &   sixteen*phi(i-1,j  ,k  ,n) - phi(i-2,j  ,k  ,n)
     & + sixteen*phi(i+1,j  ,k  ,n) - phi(i+2,j  ,k  ,n)
     & + sixteen*phi(i  ,j-1,k  ,n) - phi(i  ,j-2,k  ,n)
     & + sixteen*phi(i  ,j+1,k  ,n) - phi(i  ,j+2,k  ,n)
     & + sixteen*phi(i  ,j  ,k-1,n) - phi(i  ,j  ,k-2,n)
     & + sixteen*phi(i  ,j  ,k+1,n) - phi(i  ,j  ,k+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,j,k,n) )
     &       * twelfth * dxinv
          lofphi = alpha*phi(i,j,k,n) + beta*lap
          res(ii,jj,kk,n) = res(ii,jj,kk,n)
     &                            + (rhs(i,j,k,n) - lofphi) / denom
        
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine GSRBLAPLACIAN4(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2
     &           ,irhshi0,irhshi1,irhshi2
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,dx
     &           ,tmp
     &           ,itmplo0,itmplo1,itmplo2
     &           ,itmphi0,itmphi1,itmphi2
     &           ,ntmpcomp
     &           ,redBlack
     &           )

      implicit none
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T dx
      integer ntmpcomp
      integer itmplo0,itmplo1,itmplo2
      integer itmphi0,itmphi1,itmphi2
      REAL_T tmp(
     &           itmplo0:itmphi0,
     &           itmplo1:itmphi1,
     &           itmplo2:itmphi2,
     &           0:ntmpcomp-1)
      integer redBlack
      REAL_T dx2t, thD
      integer i,j,k
      integer n,ncomp,indtot,imin,imax,red,black
      red = 0
      black = 1
      dx2t = twelve*dx*dx
      thD  = thirtieth/CH_SPACEDIM
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         if (redBlack .eq. red) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j + k
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               tmp(i,j,k,n) = thD*( 
     &           sixteen*phi(i+1,j,k,n) - phi(i+2,j,k,n)
     &         + sixteen*phi(i-1,j,k,n) - phi(i-2,j,k,n)
     &         + sixteen*phi(i,j+1,k,n) - phi(i,j+2,k,n)
     &         + sixteen*phi(i,j-1,k,n) - phi(i,j-2,k,n)
     &         + sixteen*phi(i,j,k+1,n) - phi(i,j,k+2,n)
     &         + sixteen*phi(i,j,k-1,n) - phi(i,j,k-2,n) 
     &         - dx2t*rhs(i,j,k,n) )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else if (redBlack .eq. black) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j + k
            imin = imin + abs(mod(indtot+black, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,j,k,n) = thD*( 
     &           sixteen*tmp(i+1,j,k,n) - tmp(i+2,j,k,n)
     &         + sixteen*tmp(i-1,j,k,n) - tmp(i-2,j,k,n)
     &         + sixteen*tmp(i,j+1,k,n) - tmp(i,j+2,k,n)
     &         + sixteen*tmp(i,j-1,k,n) - tmp(i,j-2,k,n)
     &         + sixteen*tmp(i,j,k+1,n) - tmp(i,j,k+2,n)
     &         + sixteen*tmp(i,j,k-1,n) - tmp(i,j,k-2,n) 
     &         - dx2t*rhs(i,j,k,n) )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j + k
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,j,k,n) = tmp(i,j,k,n)
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else
           call MAYDAYERROR()
        end if
      enddo
      return
      end
      subroutine SORLAPLACIAN4(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2
     &           ,irhshi0,irhshi1,irhshi2
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,dx
     &           )

      implicit none
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T dx
      REAL_T dx2t, thD, tmp, omega
      integer i,j,k
      integer n,ncomp
      dx2t = twelve*dx*dx
      thD  = thirtieth/CH_SPACEDIM
      omega = 0.47
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

         tmp = thD*( 
     &        sixteen*phi(i+1,j,k,n) - phi(i+2,j,k,n)
     &        + sixteen*phi(i-1,j,k,n) - phi(i-2,j,k,n)
     &        + sixteen*phi(i,j+1,k,n) - phi(i,j+2,k,n)
     &        + sixteen*phi(i,j-1,k,n) - phi(i,j-2,k,n)
     &        + sixteen*phi(i,j,k+1,n) - phi(i,j,k+2,n)
     &        + sixteen*phi(i,j,k-1,n) - phi(i,j,k-2,n) 
     &        - dx2t*rhs(i,j,k,n) )
         phi(i,j,k,n) = omega*tmp
     &        + (one-omega)*phi(i,j,k,n)
         
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine GSRBHELMHOLTZ4(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2
     &           ,irhshi0,irhshi1,irhshi2
     &           ,nrhscomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,dx
     &           ,tmp
     &           ,itmplo0,itmplo1,itmplo2
     &           ,itmphi0,itmphi1,itmphi2
     &           ,ntmpcomp
     &           ,alpha
     &           ,beta
     &           ,redBlack
     &           )

      implicit none
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer nrhscomp
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2,
     &           0:nrhscomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T dx
      integer ntmpcomp
      integer itmplo0,itmplo1,itmplo2
      integer itmphi0,itmphi1,itmphi2
      REAL_T tmp(
     &           itmplo0:itmphi0,
     &           itmplo1:itmphi1,
     &           itmplo2:itmphi2,
     &           0:ntmpcomp-1)
      REAL_T alpha
      REAL_T beta
      integer redBlack
      REAL_T dx2t, lambda, lap, dxinv, helm
      integer i,j,k
      integer n,ncomp,indtot,imin,imax,red,black
      red = 0
      black = 1
      dx2t = twelve*dx*dx
      dxinv = one/(dx*dx)
      lambda = one/(alpha - beta*thirty*CH_SPACEDIM*twelfth*dxinv)
      lambda = lambda*(0.60)
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
         if (redBlack .eq. red) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j + k
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
          lap = ( 
     &   sixteen*phi(i-1,j  ,k  ,n) - phi(i-2,j  ,k  ,n)
     & + sixteen*phi(i+1,j  ,k  ,n) - phi(i+2,j  ,k  ,n)
     & + sixteen*phi(i  ,j-1,k  ,n) - phi(i  ,j-2,k  ,n)
     & + sixteen*phi(i  ,j+1,k  ,n) - phi(i  ,j+2,k  ,n)
     & + sixteen*phi(i  ,j  ,k-1,n) - phi(i  ,j  ,k-2,n)
     & + sixteen*phi(i  ,j  ,k+1,n) - phi(i  ,j  ,k+2,n)
     &                     -(thirty*CH_SPACEDIM)*phi(i,j,k,n) )
     &       * twelfth * dxinv
          helm = alpha*phi(i,j,k,n) + beta*lap
          tmp(i,j,k,n) = phi(i,j,k,n) +
     &      lambda*( rhs(i,j,k,n) - helm )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else if (redBlack .eq. black) then
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j + k
            imin = imin + abs(mod(indtot+black, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               lap = ( 
     &           sixteen*tmp(i+1,j,k,n) - tmp(i+2,j,k,n)
     &         + sixteen*tmp(i-1,j,k,n) - tmp(i-2,j,k,n)
     &         + sixteen*tmp(i,j+1,k,n) - tmp(i,j+2,k,n)
     &         + sixteen*tmp(i,j-1,k,n) - tmp(i,j-2,k,n)
     &         + sixteen*tmp(i,j,k+1,n) - tmp(i,j,k+2,n)
     &         + sixteen*tmp(i,j,k-1,n) - tmp(i,j,k-2,n) 
     &                     -(thirty*CH_SPACEDIM)*tmp(i,j,k,n) )
     &       * twelfth * dxinv
               helm = alpha*tmp(i,j,k,n) + beta*lap
               phi(i,j,k,n) = tmp(i,j,k,n) +
     &              lambda*( rhs(i,j,k,n) - helm )
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
#if CH_SPACEDIM > 2
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j + k
            imin = imin + abs(mod(indtot+red, 2))
            imax = iregionhi0
            do i = imin, imax, 2
               phi(i,j,k,n) = tmp(i,j,k,n)
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM > 2
        enddo
#endif
        else
           call MAYDAYERROR()
        end if
      enddo
      return
      end
      subroutine NEWGETFLUX4(
     &           flux
     &           ,ifluxlo0,ifluxlo1,ifluxlo2
     &           ,ifluxhi0,ifluxhi1,ifluxhi2
     &           ,nfluxcomp
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           ,beta_dx
     &           ,a_idir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2
      integer ifluxhi0,ifluxhi1,ifluxhi2
      REAL_T flux(
     &           ifluxlo0:ifluxhi0,
     &           ifluxlo1:ifluxhi1,
     &           ifluxlo2:ifluxhi2,
     &           0:nfluxcomp-1)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      REAL_T beta_dx
      integer a_idir
      INTEGER ncomp,n
      integer ii, jj, kk
      integer i , j , k 
      ncomp = nphicomp
      
      ii = CHF_ID(a_idir, 0)
      jj = CHF_ID(a_idir, 1)
      kk = CHF_ID(a_idir, 2)
      do n = 0, ncomp-1
          
      do k = iboxlo2,iboxhi2
      do j = iboxlo1,iboxhi1
      do i = iboxlo0,iboxhi0

          flux(i,j,k,n) = beta_dx * twelfth *
     &        ( fifteen*phi(i,j,k,n)
     &           + phi(i-2*ii,j-2*jj,k-2*kk,n)
     &           - phi(i+ii,j+jj,k+kk,n)
     &           - fifteen*phi(i-ii,j-jj,k-kk,n) )
          
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PROLONGLINEAR(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,ifineBoxlo0,ifineBoxlo1,ifineBoxlo2
     &           ,ifineBoxhi0,ifineBoxhi1,ifineBoxhi2
     &           ,icrseBoxlo0,icrseBoxlo1,icrseBoxlo2
     &           ,icrseBoxhi0,icrseBoxhi1,icrseBoxhi2
     &           ,r
     &           )

      implicit none
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer ncoarsecomp
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           icoarselo2:icoarsehi2,
     &           0:ncoarsecomp-1)
      integer ifineBoxlo0,ifineBoxlo1,ifineBoxlo2
      integer ifineBoxhi0,ifineBoxhi1,ifineBoxhi2
      integer icrseBoxlo0,icrseBoxlo1,icrseBoxlo2
      integer icrseBoxhi0,icrseBoxhi1,icrseBoxhi2
      integer r
      INTEGER ncomp, n
      integer i ,j ,k 
      integer ic,jc,kc
      ncomp = nphicomp
      do n = 0, ncomp-1
          
      do k = ifineBoxlo2,ifineBoxhi2
      do j = ifineBoxlo1,ifineBoxhi1
      do i = ifineBoxlo0,ifineBoxhi0

          
           ic = i/r
           jc = j/r
           kc = k/r
           phi(i,j,k,n) =  phi(i,j,k,n) +
     &          coarse(ic,jc,kc,n)
          
           if (ic.ne.icrseBoxhi0 .and.
     &         (ic*r.lt.i .or. ic.eq.icrseBoxlo0)) then
              phi(i,j,k,n) =  phi(i,j,k,n) +
     &             (coarse(ic+1,jc,kc,n)
     &              - coarse(ic,jc,kc,n))/r*(i+half-ic*r-half*r)
           else
              phi(i,j,k,n) =  phi(i,j,k,n) +
     &             (- coarse(ic-1,jc,kc,n)
     &              + coarse(ic,jc,kc,n))/r*(i+half-ic*r-half*r)
           endif
           if (jc.ne.icrseBoxhi1 .and.
     &         (jc*r.lt.j .or. jc.eq.icrseBoxlo1)) then
              phi(i,j,k,n) =  phi(i,j,k,n) +
     &             (coarse(ic,jc+1,kc,n)
     &              - coarse(ic,jc,kc,n))/r*(j+half-jc*r-half*r)
           else
              phi(i,j,k,n) =  phi(i,j,k,n) +
     &             (- coarse(ic,jc-1,kc,n)
     &              + coarse(ic,jc,kc,n))/r*(j+half-jc*r-half*r)
           endif
           if (kc.ne.icrseBoxhi2 .and.
     &         (kc*r.lt.k .or. kc.eq.icrseBoxlo2)) then
              phi(i,j,k,n) =  phi(i,j,k,n) +
     &             (coarse(ic,jc,kc+1,n)
     &              - coarse(ic,jc,kc,n))/r*(k+half-kc*r-half*r)
           else
              phi(i,j,k,n) =  phi(i,j,k,n) +
     &             (- coarse(ic,jc,kc-1,n)
     &              + coarse(ic,jc,kc,n))/r*(k+half-kc*r-half*r)
           endif
          
          
      enddo
      enddo
      enddo
      enddo
      return
      end
