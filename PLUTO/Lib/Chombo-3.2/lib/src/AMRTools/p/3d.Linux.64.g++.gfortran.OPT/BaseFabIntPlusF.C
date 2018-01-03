#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine BASEFABINTPLUS(
     &           sum
     &           ,isumlo0,isumlo1,isumlo2
     &           ,isumhi0,isumhi1,isumhi2
     &           ,nsumcomp
     &           ,piece
     &           ,ipiecelo0,ipiecelo1,ipiecelo2
     &           ,ipiecehi0,ipiecehi1,ipiecehi2
     &           ,npiececomp
     &           ,ibxlo0,ibxlo1,ibxlo2
     &           ,ibxhi0,ibxhi1,ibxhi2
     &           )

      implicit none
      integer nsumcomp
      integer isumlo0,isumlo1,isumlo2
      integer isumhi0,isumhi1,isumhi2
      integer sum(
     &           isumlo0:isumhi0,
     &           isumlo1:isumhi1,
     &           isumlo2:isumhi2,
     &           0:nsumcomp-1)
      integer npiececomp
      integer ipiecelo0,ipiecelo1,ipiecelo2
      integer ipiecehi0,ipiecehi1,ipiecehi2
      integer piece(
     &           ipiecelo0:ipiecehi0,
     &           ipiecelo1:ipiecehi1,
     &           ipiecelo2:ipiecehi2,
     &           0:npiececomp-1)
      integer ibxlo0,ibxlo1,ibxlo2
      integer ibxhi0,ibxhi1,ibxhi2
      integer i0,i1,i2
      integer var
      do var = 0, nsumcomp-1
         
      do i2 = ibxlo2,ibxhi2
      do i1 = ibxlo1,ibxhi1
      do i0 = ibxlo0,ibxhi0

            sum(i0,i1,i2, var) = sum(i0,i1,i2, var) +
     &        piece(i0,i1,i2, var)
         
      enddo
      enddo
      enddo
      enddo
      return
      end
