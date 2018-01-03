#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine DOTPRODUCT(
     &           dotprodout
     &           ,afab
     &           ,iafablo0,iafablo1,iafablo2
     &           ,iafabhi0,iafabhi1,iafabhi2
     &           ,nafabcomp
     &           ,bfab
     &           ,ibfablo0,ibfablo1,ibfablo2
     &           ,ibfabhi0,ibfabhi1,ibfabhi2
     &           ,nbfabcomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      REAL_T dotprodout
      integer nafabcomp
      integer iafablo0,iafablo1,iafablo2
      integer iafabhi0,iafabhi1,iafabhi2
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           iafablo2:iafabhi2,
     &           0:nafabcomp-1)
      integer nbfabcomp
      integer ibfablo0,ibfablo1,ibfablo2
      integer ibfabhi0,ibfabhi1,ibfabhi2
      REAL_T bfab(
     &           ibfablo0:ibfabhi0,
     &           ibfablo1:ibfabhi1,
     &           ibfablo2:ibfabhi2,
     &           0:nbfabcomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer startcomp
      integer endcomp
      integer i0,i1,i2
      integer nv
      dotprodout = 0
      do nv=startcomp,endcomp,1
         
      do i2 = iregionlo2,iregionhi2
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

         dotprodout = dotprodout +
     &        afab(i0,i1,i2,nv)*
     &        bfab(i0,i1,i2,nv)
         
      enddo
      enddo
      enddo
      enddo
      return
      end
