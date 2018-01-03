#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine GSMCAMRPOP(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2
     &           ,irhshi0,irhshi1,irhshi2
     &           ,icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
     &           ,icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
     &           ,dx
     &           ,alpha
     &           ,beta
     &           )

      implicit none
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2)
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2)
      integer icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
      integer icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
      REAL_T dx
      REAL_T alpha
      REAL_T beta
      real_t lambda, dxinv, sum_b, lphi
      integer i,j,k
      integer idir
      dxinv = one/(dx*dx)
      sum_b = 0.0
      do idir = 0, CH_SPACEDIM-1
         sum_b = sum_b + two*dxinv
      enddo
      lambda = one/(alpha - beta*sum_b)
      
      do k = icoloredBoxlo2,icoloredBoxhi2,2
      do j = icoloredBoxlo1,icoloredBoxhi1,2
      do i = icoloredBoxlo0,icoloredBoxhi0,2

        lphi = 
     &     (    phi(i+1,j  ,k  )
     &     +    phi(i-1,j  ,k  )
     $     -two*phi(i  ,j  ,k  )) 
     $     +(   phi(i  ,j+1,k  )
     &     +    phi(i  ,j-1,k  )
     $     -two*phi(i  ,j  ,k  )) 
     $     +(   phi(i  ,j  ,k+1)
     &     +    phi(i  ,j  ,k-1)
     $     -two*phi(i  ,j  ,k  )) 
        lphi = lphi*dxinv
        phi(i,j,k) =
     $       phi(   i,j,k) +
     &       lambda*(   rhs(   i,j,k) - lphi)
      
      enddo
      enddo
      enddo
      return
      end
      subroutine GSRBHELMHOLTZ(
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
      REAL_T alpha
      REAL_T beta
      integer redBlack
      REAL_T lambda, dxinv, sum_b, lphi, helmop
      integer i,j,k
      integer n,ncomp,idir,indtot,imin,imax
      dxinv = one/(dx*dx)
      sum_b = 0.0
      do idir = 0, CH_SPACEDIM-1
         sum_b = sum_b + two*dxinv
      enddo
      lambda = -one/(alpha - beta*sum_b)
      ncomp = nphicomp
      if (ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
#if CH_SPACEDIM==3
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j  + k 
            imin = imin + abs(mod(indtot + redBlack, 2))
            imax = iregionhi0
            do i = imin, imax, 2
              
              lphi =   (phi(i+1,j,k,n)
     &             +    phi(i-1,j,k,n)
     &             +    phi(i,j+1,k,n)
     &             +    phi(i,j-1,k,n)
     &             +    phi(i,j,k+1,n)
     &             +    phi(i,j,k-1,n)
     &             -two*CH_SPACEDIM*phi(i,j,k,n))*dxinv
              helmop = alpha*phi(i,j,k,n) + beta*lphi
              phi(i,j,k,n) = phi(i,j,k,n) +
     &             lambda*(helmop - rhs(i,j,k,n))
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM==3
        enddo
#endif
      enddo
      return
      end
      subroutine GSRBLAPLACIAN(
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
      integer redBlack
      REAL_T lambda, dxinv, sum_b, lphi, lap
      integer i,j,k
      integer n,ncomp,idir,indtot,imin,imax
      dxinv = one/(dx*dx)
      sum_b = 0.0
      do idir = 0, CH_SPACEDIM-1
         sum_b = sum_b + two*dxinv
      enddo
      lambda = -one/sum_b
      ncomp = nphicomp
      if (ncomp .ne. nrhscomp) then
         call MAYDAYERROR()
      endif
      do n = 0, ncomp - 1
#if CH_SPACEDIM==3
        do k=iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM > 1
          do j=iregionlo1, iregionhi1
#endif
            imin = iregionlo0
            indtot = imin + j  + k 
            imin = imin + abs(mod(indtot + redBlack, 2))
            imax = iregionhi0
            do i = imin, imax, 2
#if 0
              
              lphi =  ((phi(i+1,j,k,n)
     &                   - phi(i,j,k,n))
     &                   -(phi(i,j,k,n)
     &                   - phi(i-1,j,k,n)))*dxinv
     &                +  ((phi(i,j+1,k,n)
     &                   - phi(i,j,k,n))
     &                   -(phi(i,j,k,n)
     &                   - phi(i,j-1,k,n)))*dxinv
     &                +  ((phi(i,j,k+1,n)
     &                   - phi(i,j,k,n))
     &                   -(phi(i,j,k,n)
     &                   - phi(i,j,k-1,n)))*dxinv 
#else
              lphi = ( 
     &                  phi(i+1,j,k,n)
     &             +    phi(i-1,j,k,n)
     &             +    phi(i,j+1,k,n)
     &             +    phi(i,j-1,k,n)
     &             +    phi(i,j,k+1,n)
     &             +    phi(i,j,k-1,n)
     &             ) * dxinv
#endif
              phi(i,j,k,n) = lambda*(rhs(i,j,k,n)-lphi)
            enddo
#if CH_SPACEDIM > 1
          enddo
#endif
#if CH_SPACEDIM==3
        enddo
#endif
      enddo
      return
      end
      subroutine GSRBLAZY(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,lphi
     &           ,ilphilo0,ilphilo1,ilphilo2
     &           ,ilphihi0,ilphihi1,ilphihi2
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2
     &           ,irhshi0,irhshi1,irhshi2
     &           ,icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
     &           ,icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
     &           ,alpha
     &           ,beta
     &           ,dx
     &           )

      implicit none
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2)
      integer ilphilo0,ilphilo1,ilphilo2
      integer ilphihi0,ilphihi1,ilphihi2
      REAL_T lphi(
     &           ilphilo0:ilphihi0,
     &           ilphilo1:ilphihi1,
     &           ilphilo2:ilphihi2)
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2)
      integer icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
      integer icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
      REAL_T alpha
      REAL_T beta
      REAL_T dx
      integer i,j,k, idir
      real_t dxinv, sum_b, lambda
      dxinv = one/(dx*dx)
      sum_b = 0.0
      do idir = 0, CH_SPACEDIM-1
         sum_b = sum_b + two*dxinv
      enddo
      lambda = -one/(alpha - beta*sum_b)
      
      do k = icoloredBoxlo2,icoloredBoxhi2,2
      do j = icoloredBoxlo1,icoloredBoxhi1,2
      do i = icoloredBoxlo0,icoloredBoxhi0,2

      phi(i,j,k) =
     $     phi(   i,j,k) -
     &     lambda*(
     $     rhs(   i,j,k) -
     $     lphi(  i,j,k))
      
      enddo
      enddo
      enddo
      return
      end
      subroutine AMRPMULTICOLOR(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,rhs
     &           ,irhslo0,irhslo1,irhslo2
     &           ,irhshi0,irhshi1,irhshi2
     &           ,weight
     &           ,alpha
     &           ,beta
     &           ,dx
     &           ,icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
     &           ,icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
     &           )

      implicit none
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2)
      integer irhslo0,irhslo1,irhslo2
      integer irhshi0,irhshi1,irhshi2
      REAL_T rhs(
     &           irhslo0:irhshi0,
     &           irhslo1:irhshi1,
     &           irhslo2:irhshi2)
      REAL_T weight
      REAL_T alpha
      REAL_T beta
      REAL_T dx(0:2)
      integer icoloredboxlo0,icoloredboxlo1,icoloredboxlo2
      integer icoloredboxhi0,icoloredboxhi1,icoloredboxhi2
      integer i,j,k
      real_t laplphi, dx0,dx1,dx2
      dx0 = beta/(dx(0) * dx(0))
                dx1 = beta/(dx(1) * dx(1))
                dx2 = beta/(dx(2) * dx(2))
      
      do k = icoloredBoxlo2,icoloredBoxhi2,2
      do j = icoloredBoxlo1,icoloredBoxhi1,2
      do i = icoloredBoxlo0,icoloredBoxhi0,2

        laplphi = 
     &     (    phi(i+1,j  ,k  )
     &     +    phi(i-1,j  ,k  )
     $     -two*phi(i  ,j  ,k  ))*dx0 
     $     +(   phi(i  ,j+1,k  )
     &     +    phi(i  ,j-1,k  )
     $     -two*phi(i  ,j  ,k  ))*dx1 
     $     +(   phi(i  ,j  ,k+1)
     &     +    phi(i  ,j  ,k-1)
     $     -two*phi(i  ,j  ,k  ))*dx2 
        laplphi = laplphi + alpha * phi(i,j,k)
        phi(i,j,k) = phi(i,j,k) +
     &     weight*(rhs(i,j,k) - laplphi)
      
      enddo
      enddo
      enddo
      return
      end
      subroutine OPERATORLAP(
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
      if (ncomp .ne. nlofphicomp) then
         call MAYDAYERROR()
      endif
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
        
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

          lap = (  phi(i+1,j  ,k  ,n)
     &                     + phi(i-1,j  ,k  ,n) 
     &                     + phi(i  ,j+1,k  ,n)
     &                     + phi(i  ,j-1,k  ,n) 
     &                     + phi(i  ,j  ,k+1,n)
     &                     + phi(i  ,j  ,k-1,n) 
     &                     -(2*CH_SPACEDIM)*phi(i,j,k,n) )
     &       * dxinv
          lofphi(i,j,k,n) = alpha*phi(i,j,k,n)+beta*lap
        
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine OPERATORLAPRES(
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
      REAL_T dxinv, lap
      integer n,ncomp
      integer i,j,k
      ncomp = nphicomp
      dxinv = one/(dx*dx)
      do n = 0, ncomp-1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

          lap = (  phi(i+1,j  ,k  ,n)
     &                     + phi(i-1,j  ,k  ,n) 
     &                     + phi(i  ,j+1,k  ,n)
     &                     + phi(i  ,j-1,k  ,n) 
     &                     + phi(i  ,j  ,k+1,n)
     &                     + phi(i  ,j  ,k-1,n) 
     &                     -(2*CH_SPACEDIM)*phi(i,j,k,n) )
     &       * dxinv
         r(i,j,k,n) = -alpha*phi(i,j,k,n) -beta*lap +
     &       rhs(i,j,k,n)
         
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine RESTRICTRES(
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
      REAL_T denom,dxinv,lofphi
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
          lofphi = alpha * phi(i,j,k,n)
     &           + beta  *
     &              ( phi(i+1,j  ,k  ,n)
     &                        + phi(i-1,j  ,k  ,n) 
     &                        + phi(i  ,j+1,k  ,n)
     &                        + phi(i  ,j-1,k  ,n) 
     &                        + phi(i  ,j  ,k+1,n)
     &                        + phi(i  ,j  ,k-1,n) 
     &                        - phi(i  ,j  ,k  ,n) * 2 * CH_SPACEDIM
     &              ) * dxinv
          res(ii,jj,kk,n) = res(ii,jj,kk,n)
     &                            + (rhs(i,j,k,n) - lofphi) / denom
        
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PROLONG(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,m
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
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer m
      INTEGER ncomp, n
      integer i,j,k
      integer ii,jj,kk
      ncomp = nphicomp
      do n = 0, ncomp-1
          
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

          
          ii = i/m
          jj = j/m
          kk = k/m
          phi(i,j,k,n) =  phi(i,j,k,n) +
     &        coarse(ii,jj,kk,n)
         
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PROLONG_2(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,m
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
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer m
      INTEGER ncomp, n
      integer i,j,k
      integer offs(CH_SPACEDIM)
      integer ic,jc,kc
      real_t f0, den, fx(CH_SPACEDIM)
      den = one/(4**CH_SPACEDIM)
      
      fx(1) = three*den
      fx(2) = three**2*den
      fx(3) = three**3*den
      f0 = one*den
      ncomp = nphicomp
      
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

        
        ic = i/m
        jc = j/m
        kc = k/m
        
        offs(1) = 2*mod(i,2) - 1
        offs(2) = 2*mod(j,2) - 1
        offs(3) = 2*mod(k,2) - 1
        do n = 0, ncomp-1
          phi(i,j,k,n) = phi(i,j,k,n)
     $      + fx(CH_SPACEDIM)*
     $        coarse(ic,jc,kc,n)
     $      + f0*coarse(ic+offs(1),jc+offs(2),kc+offs(3),n) 
#if CH_SPACEDIM > 1
          phi(i,j,k,n) = phi(i,j,k,n) 
     $      + fx(CH_SPACEDIM-1)*
     $        (
     $         coarse(ic+offs(1),jc,kc,n)  
     $       + coarse(ic,jc+offs(2),kc,n)  
     $       + coarse(ic,jc,kc+offs(3),n)  )
#if CH_SPACEDIM > 2
          phi(i,j,k,n) = phi(i,j,k,n) 
     $      + fx(CH_SPACEDIM-2)*
     $        (
     $         coarse(ic+offs(1),jc+offs(2),kc,n)  
     $       + coarse(ic,jc+offs(2),kc+offs(3),n)  
     $       + coarse(ic+offs(1),jc,kc+offs(3),n)  )
#endif
#endif
        enddo
      
      enddo
      enddo
      enddo
      return
      end
      subroutine NEWGETFLUX(
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

          flux(i,j,k,n) =
     &        (phi(i,j,k,n)-
     &         phi(i-ii,j-jj,k-kk,n))*beta_dx
          
      enddo
      enddo
      enddo
      enddo
      return
      end
