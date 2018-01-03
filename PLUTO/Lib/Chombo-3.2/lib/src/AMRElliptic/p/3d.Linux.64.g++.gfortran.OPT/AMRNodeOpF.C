#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine NODEOPLAP(
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
      REAL_T dxinv2, lphi
      integer var, ncomp
      integer  i, j, k
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinv2 = one / (dx*dx)
      do var = 0, ncomp-1
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

            
            lphi = ( (phi(i+1 , j  , k  , var)
     &              - phi(i   , j  , k  , var) )
     &            -  (phi(i   , j  , k  , var)
     &              - phi(i-1 , j  , k  , var) ) ) * dxinv2 
     &         +   ( (phi(i   , j+1, k  , var)
     &              - phi(i   , j  , k  , var) )
     &            -  (phi(i   , j  , k  , var)
     &              - phi(i   , j-1, k  , var) ) ) * dxinv2 
     &         +   ( (phi(i   , j  , k+1, var)
     &              - phi(i   , j  , k  , var) )
     &            -  (phi(i   , j  , k  , var)
     &              - phi(i   , j  , k-1, var) ) ) * dxinv2 
            lofphi(i, j, k, var) =  lphi
         
      enddo
      enddo
      enddo
      end do
      return
      end
      subroutine NODEOPLAPPOINT(
     &           lofphi
     &           ,ilofphilo0,ilofphilo1,ilofphilo2
     &           ,ilofphihi0,ilofphihi1,ilofphihi2
     &           ,nlofphicomp
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,pt
     &           ,dx
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
      integer pt(0:2)
      REAL_T dx
      REAL_T dxinv2, lphi
      integer var, ncomp
      integer  i, j, k 
      ncomp = nphicomp
      if(ncomp .ne. nlofphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinv2 = one / (dx*dx)
      
      i = pt(0) 
      j = pt(1) 
      k = pt(2) 
      do var = 0, ncomp-1
         
         lphi = ( (phi(i+1 , j  , k  , var)
     &           - phi(i   , j  , k  , var) )
     &         -  (phi(i   , j  , k  , var)
     &           - phi(i-1 , j  , k  , var) ) ) * dxinv2 
     &      +   ( (phi(i   , j+1, k  , var)
     &           - phi(i   , j  , k  , var) )
     &         -  (phi(i   , j  , k  , var)
     &           - phi(i   , j-1, k  , var) ) ) * dxinv2 
     &      +   ( (phi(i   , j  , k+1, var)
     &           - phi(i   , j  , k  , var) )
     &         -  (phi(i   , j  , k  , var)
     &           - phi(i   , j  , k-1, var) ) ) * dxinv2 
         lofphi(i, j, k, var) =  lphi
      end do
      return
      end
      subroutine NODEGRAD(
     &           grdphi
     &           ,igrdphilo0,igrdphilo1,igrdphilo2
     &           ,igrdphihi0,igrdphihi1,igrdphihi2
     &           ,ngrdphicomp
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,dx
     &           )

      implicit none
      integer ngrdphicomp
      integer igrdphilo0,igrdphilo1,igrdphilo2
      integer igrdphihi0,igrdphihi1,igrdphihi2
      REAL_T grdphi(
     &           igrdphilo0:igrdphihi0,
     &           igrdphilo1:igrdphihi1,
     &           igrdphilo2:igrdphihi2,
     &           0:ngrdphicomp-1)
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
      REAL_T dxinvh
      integer var, ncomp, gbase
      integer  i, j, k 
      ncomp = nphicomp
      if (CH_SPACEDIM * ncomp .ne. ngrdphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinvh = half / dx
      do var = 0, ncomp-1
         gbase = CH_SPACEDIM * var
         
      do k = iregionlo2,iregionhi2
      do j = iregionlo1,iregionhi1
      do i = iregionlo0,iregionhi0

            
            grdphi(i, j, k, gbase) =
     &           ( phi(i+1 , j   , k   , var)
     &           - phi(i-1 , j   , k   , var) ) * dxinvh 
            grdphi(i, j, k, gbase + 1) =
     &           ( phi(i   , j+1 , k   , var)
     &           - phi(i   , j-1 , k   , var) ) * dxinvh 
            grdphi(i, j, k, gbase + 2) =
     &           ( phi(i   , j   , k+1 , var)
     &           - phi(i   , j   , k-1 , var) ) * dxinvh 
         
      enddo
      enddo
      enddo
      end do
      return
      end
      subroutine NODEGRADPOINT(
     &           grdphi
     &           ,igrdphilo0,igrdphilo1,igrdphilo2
     &           ,igrdphihi0,igrdphihi1,igrdphihi2
     &           ,ngrdphicomp
     &           ,phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,pt
     &           ,dx
     &           )

      implicit none
      integer ngrdphicomp
      integer igrdphilo0,igrdphilo1,igrdphilo2
      integer igrdphihi0,igrdphihi1,igrdphihi2
      REAL_T grdphi(
     &           igrdphilo0:igrdphihi0,
     &           igrdphilo1:igrdphihi1,
     &           igrdphilo2:igrdphihi2,
     &           0:ngrdphicomp-1)
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer pt(0:2)
      REAL_T dx
      REAL_T dxinvh
      integer var, ncomp, gbase
      integer  i, j, k 
      ncomp = nphicomp
      if (CH_SPACEDIM * ncomp .ne. ngrdphicomp) then
         call MAYDAY_ERROR()
      endif
      dxinvh = half / dx
      
      i = pt(0) 
      j = pt(1) 
      k = pt(2) 
      do var = 0, ncomp-1
         gbase = CH_SPACEDIM * var
            
            grdphi(i, j, k, gbase) =
     &           ( phi(i+1 , j   , k   , var)
     &           - phi(i-1 , j   , k   , var) ) * dxinvh 
            grdphi(i, j, k, gbase + 1) =
     &           ( phi(i   , j+1 , k   , var)
     &           - phi(i   , j-1 , k   , var) ) * dxinvh 
            grdphi(i, j, k, gbase + 2) =
     &           ( phi(i   , j   , k+1 , var)
     &           - phi(i   , j   , k-1 , var) ) * dxinvh 
      end do
      return
      end
      subroutine NODEGSRBLEVELLAP(
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
      REAL_T lambda
      REAL_T dxinv2, lphi
      integer i, j, k
      integer imin, imax, var, ncomp, indtot
      dxinv2 = one/(dx*dx)
      lambda = (dx*dx) / (two*CH_SPACEDIM)
      ncomp = nphicomp
      if(ncomp .ne. nrhscomp) then
         call MAYDAY_ERROR()
      endif
      do var = 0, ncomp - 1
#if CH_SPACEDIM>=3
         do k = iregionlo2, iregionhi2
#endif
#if CH_SPACEDIM>=2
            do j = iregionlo1, iregionhi1
#endif
               imin = iregionlo0
               indtot = imin  + j  + k 
               imin = imin + mod(indtot + redBlack, 2)
               imax = iregionhi0
               do i = imin, imax, 2
#ifndef NDEBUG
                  if (mod(i +j +k, 2) .ne. redBlack) then
                     print *, 'NODEGSRBLEVELLAP:  computing ',
     &                    i,  j,  k, 
     &                    ' at pass ', redBlack
                  endif
#endif
            
            lphi = ( (phi(i+1 , j  , k  , var)
     &              + phi(i-1 , j  , k  , var) ) ) * dxinv2 
     &         +   ( (phi(i   , j+1, k  , var)
     &              + phi(i   , j-1, k  , var) ) ) * dxinv2 
     &         +   ( (phi(i   , j  , k+1, var)
     &              + phi(i   , j  , k-1, var) ) ) * dxinv2 
                  phi(i,j,k, var) =
     &                 lambda * (lphi - rhs(i,j,k, var))
               enddo
#if CH_SPACEDIM>=2
            enddo
#endif
#if CH_SPACEDIM>=3
         enddo
#endif
      enddo
      return
      end
