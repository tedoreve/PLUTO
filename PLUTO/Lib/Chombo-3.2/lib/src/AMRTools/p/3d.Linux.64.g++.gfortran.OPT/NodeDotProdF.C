#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine TRAPWEIGHTS(
     &           wt
     &           ,iwtlo0,iwtlo1,iwtlo2
     &           ,iwthi0,iwthi1,iwthi2
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,dv
     &           )

      implicit none
      integer iwtlo0,iwtlo1,iwtlo2
      integer iwthi0,iwthi1,iwthi2
      REAL_T wt(
     &           iwtlo0:iwthi0,
     &           iwtlo1:iwthi1,
     &           iwtlo2:iwthi2)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T dv
      integer i0,i1,i2
      
      do i2 = iregionlo2,iregionhi2
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

         wt(i0,i1,i2) = dv
         
         if ((i0 .eq. iregionlo0) .or.
     &        (i0 .eq. iregionhi0)) then
             wt(i0,i1,i2) = wt(i0,i1,i2) * half
         endif 
         if ((i1 .eq. iregionlo1) .or.
     &        (i1 .eq. iregionhi1)) then
             wt(i0,i1,i2) = wt(i0,i1,i2) * half
         endif 
         if ((i2 .eq. iregionlo2) .or.
     &        (i2 .eq. iregionhi2)) then
             wt(i0,i1,i2) = wt(i0,i1,i2) * half
         endif 
      
      enddo
      enddo
      enddo
      return
      end
      subroutine NODEDOTPRODUCT(
     &           dotprodout
     &           ,afab
     &           ,iafablo0,iafablo1,iafablo2
     &           ,iafabhi0,iafabhi1,iafabhi2
     &           ,nafabcomp
     &           ,bfab
     &           ,ibfablo0,ibfablo1,ibfablo2
     &           ,ibfabhi0,ibfabhi1,ibfabhi2
     &           ,nbfabcomp
     &           ,wfab
     &           ,iwfablo0,iwfablo1,iwfablo2
     &           ,iwfabhi0,iwfabhi1,iwfabhi2
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
      integer iwfablo0,iwfablo1,iwfablo2
      integer iwfabhi0,iwfabhi1,iwfabhi2
      REAL_T wfab(
     &           iwfablo0:iwfabhi0,
     &           iwfablo1:iwfabhi1,
     &           iwfablo2:iwfabhi2)
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
     &        bfab(i0,i1,i2,nv)*
     &        wfab(i0,i1,i2)
         
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine NODENORM(
     &           normout
     &           ,p
     &           ,afab
     &           ,iafablo0,iafablo1,iafablo2
     &           ,iafabhi0,iafabhi1,iafabhi2
     &           ,nafabcomp
     &           ,wfab
     &           ,iwfablo0,iwfablo1,iwfablo2
     &           ,iwfabhi0,iwfabhi1,iwfabhi2
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      REAL_T normout
      integer p
      integer nafabcomp
      integer iafablo0,iafablo1,iafablo2
      integer iafabhi0,iafabhi1,iafabhi2
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           iafablo2:iafabhi2,
     &           0:nafabcomp-1)
      integer iwfablo0,iwfablo1,iwfablo2
      integer iwfabhi0,iwfabhi1,iwfabhi2
      REAL_T wfab(
     &           iwfablo0:iwfabhi0,
     &           iwfablo1:iwfabhi1,
     &           iwfablo2:iwfabhi2)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer startcomp
      integer endcomp
      integer i0,i1,i2
      integer nv
      REAL_T pwrinv
      normout = 0
      if (p .eq. 1) then
         do nv = startcomp, endcomp, 1
            
      do i2 = iregionlo2,iregionhi2
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

               normout = normout +
     &           wfab(i0,i1,i2) *
     &           abs(afab(i0,i1,i2, nv))
            
      enddo
      enddo
      enddo
         enddo
      elseif (p .eq. 2) then
         do nv = startcomp, endcomp, 1
            
      do i2 = iregionlo2,iregionhi2
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

            normout = normout +
     &           wfab(i0,i1,i2) *
     &           afab(i0,i1,i2, nv) *
     &           afab(i0,i1,i2, nv)
            
      enddo
      enddo
      enddo
         enddo
         normout = sqrt(normout)
      else
         do nv = startcomp, endcomp, 1
            
      do i2 = iregionlo2,iregionhi2
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

               normout = normout +
     &           wfab(i0,i1,i2) *
     &           (afab(i0,i1,i2, nv)**p)
            
      enddo
      enddo
      enddo
         enddo
         pwrinv = one / p
         normout = normout**pwrinv
      endif
      return
      end
      subroutine NODEINTEGRAL(
     &           ans
     &           ,afab
     &           ,iafablo0,iafablo1,iafablo2
     &           ,iafabhi0,iafabhi1,iafabhi2
     &           ,nafabcomp
     &           ,wfab
     &           ,iwfablo0,iwfablo1,iwfablo2
     &           ,iwfabhi0,iwfabhi1,iwfabhi2
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      REAL_T ans
      integer nafabcomp
      integer iafablo0,iafablo1,iafablo2
      integer iafabhi0,iafabhi1,iafabhi2
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           iafablo2:iafabhi2,
     &           0:nafabcomp-1)
      integer iwfablo0,iwfablo1,iwfablo2
      integer iwfabhi0,iwfabhi1,iwfabhi2
      REAL_T wfab(
     &           iwfablo0:iwfabhi0,
     &           iwfablo1:iwfabhi1,
     &           iwfablo2:iwfabhi2)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer startcomp
      integer endcomp
      integer i0,i1,i2
      integer nv
      ans = 0
      do nv = startcomp, endcomp, 1
         
      do i2 = iregionlo2,iregionhi2
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

            ans = ans +
     &        wfab(i0,i1,i2) * afab(i0,i1,i2, nv)
         
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine NODEMAXNORM(
     &           normout
     &           ,afab
     &           ,iafablo0,iafablo1,iafablo2
     &           ,iafabhi0,iafabhi1,iafabhi2
     &           ,nafabcomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      REAL_T normout
      integer nafabcomp
      integer iafablo0,iafablo1,iafablo2
      integer iafabhi0,iafabhi1,iafabhi2
      REAL_T afab(
     &           iafablo0:iafabhi0,
     &           iafablo1:iafabhi1,
     &           iafablo2:iafabhi2,
     &           0:nafabcomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer startcomp
      integer endcomp
      integer i0,i1,i2
      integer nv
      REAL_T this
      normout = 0
      do nv = startcomp, endcomp, 1
         
      do i2 = iregionlo2,iregionhi2
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

            this = abs(afab(i0,i1,i2, nv))
            if (this .gt. normout) then
               normout = this
            endif
         
      enddo
      enddo
      enddo
      enddo
      return
      end
