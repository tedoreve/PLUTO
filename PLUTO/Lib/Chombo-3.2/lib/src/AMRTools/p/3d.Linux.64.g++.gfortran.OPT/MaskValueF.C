#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine MASKVALUE(
     &           mask
     &           ,imasklo0,imasklo1,imasklo2
     &           ,imaskhi0,imaskhi1,imaskhi2
     &           ,test
     &           ,itestlo0,itestlo1,itestlo2
     &           ,itesthi0,itesthi1,itesthi2
     &           ,ibxlo0,ibxlo1,ibxlo2
     &           ,ibxhi0,ibxhi1,ibxhi2
     &           ,val
     &           ,onoff
     &           )

      implicit none
      integer imasklo0,imasklo1,imasklo2
      integer imaskhi0,imaskhi1,imaskhi2
      REAL_T mask(
     &           imasklo0:imaskhi0,
     &           imasklo1:imaskhi1,
     &           imasklo2:imaskhi2)
      integer itestlo0,itestlo1,itestlo2
      integer itesthi0,itesthi1,itesthi2
      integer test(
     &           itestlo0:itesthi0,
     &           itestlo1:itesthi1,
     &           itestlo2:itesthi2)
      integer ibxlo0,ibxlo1,ibxlo2
      integer ibxhi0,ibxhi1,ibxhi2
      integer val
      integer onoff
      integer i0,i1,i2
      if (onoff .eq. 1) then
         
      do i2 = ibxlo2,ibxhi2
      do i1 = ibxlo1,ibxhi1
      do i0 = ibxlo0,ibxhi0

         if (test(i0,i1,i2) .eq. val) then
            mask(i0,i1,i2) = one
         else
            mask(i0,i1,i2) = 0
         endif
         
      enddo
      enddo
      enddo
      else
         
      do i2 = ibxlo2,ibxhi2
      do i1 = ibxlo1,ibxhi1
      do i0 = ibxlo0,ibxhi0

         if (test(i0,i1,i2) .eq. val) then
            mask(i0,i1,i2) = 0
         else
            mask(i0,i1,i2) = one
         endif
         
      enddo
      enddo
      enddo
      endif
      return
      end
