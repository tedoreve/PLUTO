#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine AVERAGE(
     &           coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,nfinecomp
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           ,refRatio
     &           ,ibreflo0,ibreflo1,ibreflo2
     &           ,ibrefhi0,ibrefhi1,ibrefhi2
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
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer refRatio
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer var
      integer ic0,ic1,ic2
      integer ip0,ip1,ip2
      integer ii0,ii1,ii2
      real_t refScale,coarseSum
      refScale = one / (refRatio**CH_SPACEDIM)
#if (CH_SPACEDIM < 4)
      if (refRatio .eq. 2) then
        do var = 0, ncoarsecomp - 1
          
      do ic2 = iboxlo2,iboxhi2
      do ic1 = iboxlo1,iboxhi1
      do ic0 = iboxlo0,iboxhi0

            
            ip0 = 2*ic0
            ip1 = 2*ic1
            ip2 = 2*ic2
            coarse(ic0,ic1,ic2,var) = refScale *
     &        (
     &          fine(ip0,ip1,ip2,var)
     &        + fine(ip0+1,ip1  ,ip2  ,var)
     &        + fine(ip0  ,ip1+1,ip2  ,var)
     &        + fine(ip0+1,ip1+1,ip2  ,var)
     &        + fine(ip0  ,ip1  ,ip2+1,var)
     &        + fine(ip0+1,ip1  ,ip2+1,var)
     &        + fine(ip0  ,ip1+1,ip2+1,var)
     &        + fine(ip0+1,ip1+1,ip2+1,var) )
          
      enddo
      enddo
      enddo
        enddo
      else if (refRatio .eq. 4) then
        do var = 0, ncoarsecomp - 1
          
      do ic2 = iboxlo2,iboxhi2
      do ic1 = iboxlo1,iboxhi1
      do ic0 = iboxlo0,iboxhi0

            
            ip0 = 4*ic0
            ip1 = 4*ic1
            ip2 = 4*ic2
            coarse(ic0,ic1,ic2,var) = refScale *
     &        (
     &           fine(ip0  ,ip1  ,ip2  ,var)
     &         + fine(ip0+1,ip1  ,ip2  ,var)
     &         + fine(ip0+2,ip1  ,ip2  ,var)
     &         + fine(ip0+3,ip1  ,ip2  ,var) 
     &         + fine(ip0  ,ip1+1,ip2  ,var)
     &         + fine(ip0+1,ip1+1,ip2  ,var)
     &         + fine(ip0+2,ip1+1,ip2  ,var)
     &         + fine(ip0+3,ip1+1,ip2  ,var)
     &         + fine(ip0  ,ip1+2,ip2  ,var)
     &         + fine(ip0+1,ip1+2,ip2  ,var)
     &         + fine(ip0+2,ip1+2,ip2  ,var)
     &         + fine(ip0+3,ip1+2,ip2  ,var)
     &         + fine(ip0  ,ip1+3,ip2  ,var)
     &         + fine(ip0+1,ip1+3,ip2  ,var)
     &         + fine(ip0+2,ip1+3,ip2  ,var)
     &         + fine(ip0+3,ip1+3,ip2  ,var) 
     &         + fine(ip0  ,ip1  ,ip2+1,var)
     &         + fine(ip0+1,ip1  ,ip2+1,var)
     &         + fine(ip0+2,ip1  ,ip2+1,var)
     &         + fine(ip0+3,ip1  ,ip2+1,var)
     &         + fine(ip0  ,ip1+1,ip2+1,var)
     &         + fine(ip0+1,ip1+1,ip2+1,var)
     &         + fine(ip0+2,ip1+1,ip2+1,var)
     &         + fine(ip0+3,ip1+1,ip2+1,var)
     &         + fine(ip0  ,ip1+2,ip2+1,var)
     &         + fine(ip0+1,ip1+2,ip2+1,var)
     &         + fine(ip0+2,ip1+2,ip2+1,var)
     &         + fine(ip0+3,ip1+2,ip2+1,var)
     &         + fine(ip0  ,ip1+3,ip2+1,var)
     &         + fine(ip0+1,ip1+3,ip2+1,var)
     &         + fine(ip0+2,ip1+3,ip2+1,var)
     &         + fine(ip0+3,ip1+3,ip2+1,var)
     &         + fine(ip0  ,ip1  ,ip2+2,var)
     &         + fine(ip0+1,ip1  ,ip2+2,var)
     &         + fine(ip0+2,ip1  ,ip2+2,var)
     &         + fine(ip0+3,ip1  ,ip2+2,var)
     &         + fine(ip0  ,ip1+1,ip2+2,var)
     &         + fine(ip0+1,ip1+1,ip2+2,var)
     &         + fine(ip0+2,ip1+1,ip2+2,var)
     &         + fine(ip0+3,ip1+1,ip2+2,var)
     &         + fine(ip0  ,ip1+2,ip2+2,var)
     &         + fine(ip0+1,ip1+2,ip2+2,var)
     &         + fine(ip0+2,ip1+2,ip2+2,var)
     &         + fine(ip0+3,ip1+2,ip2+2,var)
     &         + fine(ip0  ,ip1+3,ip2+2,var)
     &         + fine(ip0+1,ip1+3,ip2+2,var)
     &         + fine(ip0+2,ip1+3,ip2+2,var)
     &         + fine(ip0+3,ip1+3,ip2+2,var)
     &         + fine(ip0  ,ip1  ,ip2+3,var)
     &         + fine(ip0+1,ip1  ,ip2+3,var)
     &         + fine(ip0+2,ip1  ,ip2+3,var)
     &         + fine(ip0+3,ip1  ,ip2+3,var)
     &         + fine(ip0  ,ip1+1,ip2+3,var)
     &         + fine(ip0+1,ip1+1,ip2+3,var)
     &         + fine(ip0+2,ip1+1,ip2+3,var)
     &         + fine(ip0+3,ip1+1,ip2+3,var)
     &         + fine(ip0  ,ip1+2,ip2+3,var)
     &         + fine(ip0+1,ip1+2,ip2+3,var)
     &         + fine(ip0+2,ip1+2,ip2+3,var)
     &         + fine(ip0+3,ip1+2,ip2+3,var)
     &         + fine(ip0  ,ip1+3,ip2+3,var)
     &         + fine(ip0+1,ip1+3,ip2+3,var)
     &         + fine(ip0+2,ip1+3,ip2+3,var)
     &         + fine(ip0+3,ip1+3,ip2+3,var) )
          
      enddo
      enddo
      enddo
        enddo
      else
#endif
        do var = 0, ncoarsecomp - 1
          
      do ic2 = iboxlo2,iboxhi2
      do ic1 = iboxlo1,iboxhi1
      do ic0 = iboxlo0,iboxhi0

            
            ip0 = ic0*refRatio
            ip1 = ic1*refRatio
            ip2 = ic2*refRatio
            coarseSum = 0
            
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

              coarseSum = coarseSum + fine( ip0+ii0,ip1+ii1,ip2+ii2,var)
            
      enddo
      enddo
      enddo
            coarse(ic0,ic1,ic2,var) = coarseSum * refScale
          
      enddo
      enddo
      enddo
       enddo
#if (CH_SPACEDIM < 4)
      endif
#endif
      return
      end
      subroutine AVERAGEHARMONIC(
     &           coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,nfinecomp
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           ,refRatio
     &           ,ibreflo0,ibreflo1,ibreflo2
     &           ,ibrefhi0,ibrefhi1,ibrefhi2
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
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer refRatio
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer var
      integer ic0,ic1,ic2
      integer ip0,ip1,ip2
      integer ii0,ii1,ii2
      real_t refScale,coarseSum
      refScale = one / (refRatio**CH_SPACEDIM)
#if (CH_SPACEDIM < 4)
      if (refRatio .eq. 2) then
        do var = 0, ncoarsecomp - 1
          
      do ic2 = iboxlo2,iboxhi2
      do ic1 = iboxlo1,iboxhi1
      do ic0 = iboxlo0,iboxhi0

            
            ip0 = 2*ic0
            ip1 = 2*ic1
            ip2 = 2*ic2
            coarse(ic0,ic1,ic2,var) = refScale *
     &        (
     &          one/fine(ip0  ,ip1  ,ip2  ,var)
     &        + one/fine(ip0+1,ip1  ,ip2  ,var)
     &        + one/fine(ip0  ,ip1+1,ip2  ,var)
     &        + one/fine(ip0+1,ip1+1,ip2  ,var)
     &        + one/fine(ip0  ,ip1  ,ip2+1,var)
     &        + one/fine(ip0+1,ip1  ,ip2+1,var)
     &        + one/fine(ip0  ,ip1+1,ip2+1,var)
     &        + one/fine(ip0+1,ip1+1,ip2+1,var) )
            coarse(ic0,ic1,ic2,var) = one/coarse(ic0,ic1,ic2,var)
          
      enddo
      enddo
      enddo
        enddo
      else if (refRatio .eq. 4) then
        do var = 0, ncoarsecomp - 1
          
      do ic2 = iboxlo2,iboxhi2
      do ic1 = iboxlo1,iboxhi1
      do ic0 = iboxlo0,iboxhi0

            
            ip0 = 4*ic0
            ip1 = 4*ic1
            ip2 = 4*ic2
            coarse(ic0,ic1,ic2,var) = refScale *
     &        (
     &           one/fine(ip0  ,ip1  ,ip2  ,var)
     &         + one/fine(ip0+1,ip1  ,ip2  ,var)
     &         + one/fine(ip0+2,ip1  ,ip2  ,var)
     &         + one/fine(ip0+3,ip1  ,ip2  ,var) 
     &         + one/fine(ip0  ,ip1+1,ip2  ,var)
     &         + one/fine(ip0+1,ip1+1,ip2  ,var)
     &         + one/fine(ip0+2,ip1+1,ip2  ,var)
     &         + one/fine(ip0+3,ip1+1,ip2  ,var)
     &         + one/fine(ip0  ,ip1+2,ip2  ,var)
     &         + one/fine(ip0+1,ip1+2,ip2  ,var)
     &         + one/fine(ip0+2,ip1+2,ip2  ,var)
     &         + one/fine(ip0+3,ip1+2,ip2  ,var)
     &         + one/fine(ip0  ,ip1+3,ip2  ,var)
     &         + one/fine(ip0+1,ip1+3,ip2  ,var)
     &         + one/fine(ip0+2,ip1+3,ip2  ,var)
     &         + one/fine(ip0+3,ip1+3,ip2  ,var) 
     &         + one/fine(ip0  ,ip1  ,ip2+1,var)
     &         + one/fine(ip0+1,ip1  ,ip2+1,var)
     &         + one/fine(ip0+2,ip1  ,ip2+1,var)
     &         + one/fine(ip0+3,ip1  ,ip2+1,var)
     &         + one/fine(ip0  ,ip1+1,ip2+1,var)
     &         + one/fine(ip0+1,ip1+1,ip2+1,var)
     &         + one/fine(ip0+2,ip1+1,ip2+1,var)
     &         + one/fine(ip0+3,ip1+1,ip2+1,var)
     &         + one/fine(ip0  ,ip1+2,ip2+1,var)
     &         + one/fine(ip0+1,ip1+2,ip2+1,var)
     &         + one/fine(ip0+2,ip1+2,ip2+1,var)
     &         + one/fine(ip0+3,ip1+2,ip2+1,var)
     &         + one/fine(ip0  ,ip1+3,ip2+1,var)
     &         + one/fine(ip0+1,ip1+3,ip2+1,var)
     &         + one/fine(ip0+2,ip1+3,ip2+1,var)
     &         + one/fine(ip0+3,ip1+3,ip2+1,var)
     &         + one/fine(ip0  ,ip1  ,ip2+2,var)
     &         + one/fine(ip0+1,ip1  ,ip2+2,var)
     &         + one/fine(ip0+2,ip1  ,ip2+2,var)
     &         + one/fine(ip0+3,ip1  ,ip2+2,var)
     &         + one/fine(ip0  ,ip1+1,ip2+2,var)
     &         + one/fine(ip0+1,ip1+1,ip2+2,var)
     &         + one/fine(ip0+2,ip1+1,ip2+2,var)
     &         + one/fine(ip0+3,ip1+1,ip2+2,var)
     &         + one/fine(ip0  ,ip1+2,ip2+2,var)
     &         + one/fine(ip0+1,ip1+2,ip2+2,var)
     &         + one/fine(ip0+2,ip1+2,ip2+2,var)
     &         + one/fine(ip0+3,ip1+2,ip2+2,var)
     &         + one/fine(ip0  ,ip1+3,ip2+2,var)
     &         + one/fine(ip0+1,ip1+3,ip2+2,var)
     &         + one/fine(ip0+2,ip1+3,ip2+2,var)
     &         + one/fine(ip0+3,ip1+3,ip2+2,var)
     &         + one/fine(ip0  ,ip1  ,ip2+3,var)
     &         + one/fine(ip0+1,ip1  ,ip2+3,var)
     &         + one/fine(ip0+2,ip1  ,ip2+3,var)
     &         + one/fine(ip0+3,ip1  ,ip2+3,var)
     &         + one/fine(ip0  ,ip1+1,ip2+3,var)
     &         + one/fine(ip0+1,ip1+1,ip2+3,var)
     &         + one/fine(ip0+2,ip1+1,ip2+3,var)
     &         + one/fine(ip0+3,ip1+1,ip2+3,var)
     &         + one/fine(ip0  ,ip1+2,ip2+3,var)
     &         + one/fine(ip0+1,ip1+2,ip2+3,var)
     &         + one/fine(ip0+2,ip1+2,ip2+3,var)
     &         + one/fine(ip0+3,ip1+2,ip2+3,var)
     &         + one/fine(ip0  ,ip1+3,ip2+3,var)
     &         + one/fine(ip0+1,ip1+3,ip2+3,var)
     &         + one/fine(ip0+2,ip1+3,ip2+3,var)
     &         + one/fine(ip0+3,ip1+3,ip2+3,var) )
            coarse(ic0,ic1,ic2,var) = one/coarse(ic0,ic1,ic2,var)
          
      enddo
      enddo
      enddo
        enddo
      else
#endif
        do var = 0, ncoarsecomp - 1
          
      do ic2 = iboxlo2,iboxhi2
      do ic1 = iboxlo1,iboxhi1
      do ic0 = iboxlo0,iboxhi0

            
            ip0 = ic0*refRatio
            ip1 = ic1*refRatio
            ip2 = ic2*refRatio
            coarseSum = 0
            
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

              coarseSum = coarseSum + one/fine( ip0+ii0,ip1+ii1,ip2+ii2,var)
            
      enddo
      enddo
      enddo
            coarse(ic0,ic1,ic2,var) = one/(coarseSum * refScale)
          
      enddo
      enddo
      enddo
       enddo
#if (CH_SPACEDIM < 4)
      endif
#endif
      return
      end
      subroutine AVERAGEINTVECTREF(
     &           coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,nfinecomp
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           ,ref
     &           ,ibreflo0,ibreflo1,ibreflo2
     &           ,ibrefhi0,ibrefhi1,ibrefhi2
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
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer ref(0:2)
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer var
      integer ic0,ic1,ic2
      integer ip0,ip1,ip2
      integer ii0,ii1,ii2
      real_t refScale,coarseSum
      refScale = one / (
     &  ref(0)*ref(1)*ref(2))
      do var = 0, ncoarsecomp - 1
        
      do ic2 = iboxlo2,iboxhi2
      do ic1 = iboxlo1,iboxhi1
      do ic0 = iboxlo0,iboxhi0

          
          ip0 = ic0*ref(0)
          ip1 = ic1*ref(1)
          ip2 = ic2*ref(2)
          coarseSum = 0
          
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

            coarseSum = coarseSum + fine( ip0+ii0,ip1+ii1,ip2+ii2,var)
          
      enddo
      enddo
      enddo
          coarse(ic0,ic1,ic2,var) = coarseSum * refScale
        
      enddo
      enddo
      enddo
      enddo
      return
      end
