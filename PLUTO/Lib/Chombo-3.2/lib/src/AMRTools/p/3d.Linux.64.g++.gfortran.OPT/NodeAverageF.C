#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine NODEAVERAGE(
     &           coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,nfinecomp
     &           ,iblo0,iblo1,iblo2
     &           ,ibhi0,ibhi1,ibhi2
     &           ,ref_ratio
     &           ,weight
     &           ,iweightlo0,iweightlo1,iweightlo2
     &           ,iweighthi0,iweighthi1,iweighthi2
     &           ,nweightcomp
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
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer ref_ratio
      integer nweightcomp
      integer iweightlo0,iweightlo1,iweightlo2
      integer iweighthi0,iweighthi1,iweighthi2
      REAL_T weight(
     &           iweightlo0:iweighthi0,
     &           iweightlo1:iweighthi1,
     &           iweightlo2:iweighthi2,
     &           0:nweightcomp-1)
      integer var
      integer icrse0,icrse1,icrse2
      integer ifine0,ifine1,ifine2
      integer ii0,ii1,ii2
      REAL_T csum
      do var = 0, ncoarsecomp - 1
         
      do icrse2 = iblo2,ibhi2
      do icrse1 = iblo1,ibhi1
      do icrse0 = iblo0,ibhi0

            csum = 0
            
      do ii2 = iweightlo2,iweighthi2
      do ii1 = iweightlo1,iweighthi1
      do ii0 = iweightlo0,iweighthi0

               
               ifine0 = icrse0*ref_ratio + ii0 
               ifine1 = icrse1*ref_ratio + ii1 
               ifine2 = icrse2*ref_ratio + ii2 
               csum = csum + weight( ii0,ii1,ii2, 0) *
     &                 fine( ifine0,ifine1,ifine2, var )
            
      enddo
      enddo
      enddo
            coarse( icrse0,icrse1,icrse2, var ) = csum
         
      enddo
      enddo
      enddo
      end do
      return
      end
      subroutine NODEAVERAGEPOINT(
     &           coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,nfinecomp
     &           ,pcrse
     &           ,ref_ratio
     &           ,weight
     &           ,iweightlo0,iweightlo1,iweightlo2
     &           ,iweighthi0,iweighthi1,iweighthi2
     &           ,nweightcomp
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
      integer pcrse(0:2)
      integer ref_ratio
      integer nweightcomp
      integer iweightlo0,iweightlo1,iweightlo2
      integer iweighthi0,iweighthi1,iweighthi2
      REAL_T weight(
     &           iweightlo0:iweighthi0,
     &           iweightlo1:iweighthi1,
     &           iweightlo2:iweighthi2,
     &           0:nweightcomp-1)
      integer var
      integer ifine0,ifine1,ifine2
      integer ii0,ii1,ii2
      REAL_T csum, weightpt, finept
      do var = 0, ncoarsecomp - 1
         csum = 0
         
      do ii2 = iweightlo2,iweighthi2
      do ii1 = iweightlo1,iweighthi1
      do ii0 = iweightlo0,iweighthi0

            
            ifine0 = pcrse(0)*ref_ratio + ii0 
            ifine1 = pcrse(1)*ref_ratio + ii1 
            ifine2 = pcrse(2)*ref_ratio + ii2 
            weightpt = weight( ii0,ii1,ii2, 0)
            finept = fine( ifine0,ifine1,ifine2, var )
            csum = csum + weightpt*finept
         
      enddo
      enddo
      enddo
         coarse(pcrse(0),pcrse(1),pcrse(2), var ) = csum
      end do
      return
      end
      subroutine NODEAVERAGE_GETWEIGHTS(
     &           weight
     &           ,iweightlo0,iweightlo1,iweightlo2
     &           ,iweighthi0,iweighthi1,iweighthi2
     &           ,nweightcomp
     &           ,ref_ratio
     &           )

      implicit none
      integer nweightcomp
      integer iweightlo0,iweightlo1,iweightlo2
      integer iweighthi0,iweighthi1,iweighthi2
      REAL_T weight(
     &           iweightlo0:iweighthi0,
     &           iweightlo1:iweighthi1,
     &           iweightlo2:iweighthi2,
     &           0:nweightcomp-1)
      integer ref_ratio
      integer ext, nxtrm
      integer ii0,ii1,ii2
      REAL_T ref_scale
      ext = ref_ratio / 2
      ref_scale = one / (ref_ratio**CH_SPACEDIM)
      
      do ii2 = iweightlo2,iweighthi2
      do ii1 = iweightlo1,iweighthi1
      do ii0 = iweightlo0,iweighthi0

         nxtrm = 0
         
         if (iabs(ii0) .eq. ext) nxtrm = nxtrm + 1 
         if (iabs(ii1) .eq. ext) nxtrm = nxtrm + 1 
         if (iabs(ii2) .eq. ext) nxtrm = nxtrm + 1 
         weight( ii0,ii1,ii2, 0) = ref_scale * half**nxtrm
      
      enddo
      enddo
      enddo
      return
      end
