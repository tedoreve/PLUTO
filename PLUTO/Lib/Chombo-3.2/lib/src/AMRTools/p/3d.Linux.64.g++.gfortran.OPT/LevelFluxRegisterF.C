#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine INCREMENTFINE(
     &           fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,nfinecomp
     &           ,cFine
     &           ,icFinelo0,icFinelo1,icFinelo2
     &           ,icFinehi0,icFinehi1,icFinehi2
     &           ,ncFinecomp
     &           ,ifineBoxlo0,ifineBoxlo1,ifineBoxlo2
     &           ,ifineBoxhi0,ifineBoxhi1,ifineBoxhi2
     &           ,nRef
     &           ,scale
     &           ,srcStart
     &           ,destStart
     &           ,ncomp
     &           )

      implicit none
      integer nfinecomp
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           ifinelo2:ifinehi2,
     &           0:nfinecomp-1)
      integer ncFinecomp
      integer icFinelo0,icFinelo1,icFinelo2
      integer icFinehi0,icFinehi1,icFinehi2
      REAL_T cFine(
     &           icFinelo0:icFinehi0,
     &           icFinelo1:icFinehi1,
     &           icFinelo2:icFinehi2,
     &           0:ncFinecomp-1)
      integer ifineBoxlo0,ifineBoxlo1,ifineBoxlo2
      integer ifineBoxhi0,ifineBoxhi1,ifineBoxhi2
      integer nRef(0:2)
      REAL_T scale
      integer srcStart
      integer destStart
      integer ncomp
      integer i0,i1,i2
      integer ii0,ii1,ii2
      integer var, srcComp, destComp
      do var=0, ncomp-1
         srcComp = srcStart + var
         destComp = destStart + var
         
      do i2 = ifineBoxlo2,ifineBoxhi2
      do i1 = ifineBoxlo1,ifineBoxhi1
      do i0 = ifineBoxlo0,ifineBoxhi0

            
            ii0=i0/nRef(0)
            ii1=i1/nRef(1)
            ii2=i2/nRef(2)
            cFine(ii0,ii1,ii2,destComp) =
     &           cFine(ii0,ii1,ii2, destComp) +
     &           scale * fine(i0,i1,i2, srcComp)
         
      enddo
      enddo
      enddo
      enddo
      return
      end
