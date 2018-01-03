#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine AVERAGEFACE(
     &           coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,nfinecomp
     &           ,icrseBoxlo0,icrseBoxlo1,icrseBoxlo2
     &           ,icrseBoxhi0,icrseBoxhi1,icrseBoxhi2
     &           ,dir
     &           ,nRef
     &           ,refFactor
     &           ,irefBoxlo0,irefBoxlo1,irefBoxlo2
     &           ,irefBoxhi0,irefBoxhi1,irefBoxhi2
     &           )

      implicit none
      integer ncoarsecomp
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           icoarselo2:icoarsehi2,
     &           0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           ifinelo2:ifinehi2,
     &           0:nfinecomp-1)
      integer icrseBoxlo0,icrseBoxlo1,icrseBoxlo2
      integer icrseBoxhi0,icrseBoxhi1,icrseBoxhi2
      integer dir
      integer nRef
      integer refFactor
      integer irefBoxlo0,irefBoxlo1,irefBoxlo2
      integer irefBoxhi0,irefBoxhi1,irefBoxhi2
      integer ic0,ic1,ic2
      integer ifine0,ifine1,ifine2
      integer var
      integer ii0,ii1,ii2
      REAL_T crseSum, ref_scale
      ref_scale = (one/refFactor)**(CH_SPACEDIM-1)
      do var=0, ncoarsecomp-1
         
      do ic2 = icrseBoxlo2,icrseBoxhi2
      do ic1 = icrseBoxlo1,icrseBoxhi1
      do ic0 = icrseBoxlo0,icrseBoxhi0

         crseSum = 0
         
      do ii2 = irefBoxlo2,irefBoxhi2
      do ii1 = irefBoxlo1,irefBoxhi1
      do ii0 = irefBoxlo0,irefBoxhi0

         
         ifine0=ic0*nRef+ii0
         ifine1=ic1*nRef+ii1
         ifine2=ic2*nRef+ii2
            crseSum = crseSum + fine(ifine0,ifine1,ifine2,var)
            
      enddo
      enddo
      enddo
            coarse(ic0,ic1,ic2,var) = ref_scale*crseSum
          
      enddo
      enddo
      enddo
       enddo
       return
       end
      subroutine AVERAGEFACEHARMONIC(
     &           coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,nfinecomp
     &           ,icrseBoxlo0,icrseBoxlo1,icrseBoxlo2
     &           ,icrseBoxhi0,icrseBoxhi1,icrseBoxhi2
     &           ,dir
     &           ,nRef
     &           ,refFactor
     &           ,irefBoxlo0,irefBoxlo1,irefBoxlo2
     &           ,irefBoxhi0,irefBoxhi1,irefBoxhi2
     &           )

      implicit none
      integer ncoarsecomp
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           icoarselo2:icoarsehi2,
     &           0:ncoarsecomp-1)
      integer nfinecomp
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           ifinelo2:ifinehi2,
     &           0:nfinecomp-1)
      integer icrseBoxlo0,icrseBoxlo1,icrseBoxlo2
      integer icrseBoxhi0,icrseBoxhi1,icrseBoxhi2
      integer dir
      integer nRef
      integer refFactor
      integer irefBoxlo0,irefBoxlo1,irefBoxlo2
      integer irefBoxhi0,irefBoxhi1,irefBoxhi2
      integer ic0,ic1,ic2
      integer ifine0,ifine1,ifine2
      integer var
      integer ii0,ii1,ii2
      REAL_T crseSum, ref_scale
      ref_scale = (one/refFactor)**(CH_SPACEDIM-1)
      do var=0, ncoarsecomp-1
         
      do ic2 = icrseBoxlo2,icrseBoxhi2
      do ic1 = icrseBoxlo1,icrseBoxhi1
      do ic0 = icrseBoxlo0,icrseBoxhi0

         crseSum = 0
         
      do ii2 = irefBoxlo2,irefBoxhi2
      do ii1 = irefBoxlo1,irefBoxhi1
      do ii0 = irefBoxlo0,irefBoxhi0

         
         ifine0=ic0*nRef+ii0
         ifine1=ic1*nRef+ii1
         ifine2=ic2*nRef+ii2
            crseSum = crseSum + one/fine(ifine0,ifine1,ifine2,var)
            
      enddo
      enddo
      enddo
            coarse(ic0,ic1,ic2,var) = one/(ref_scale*crseSum)
          
      enddo
      enddo
      enddo
       enddo
       return
       end
