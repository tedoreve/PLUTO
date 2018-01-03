      subroutine NODEINTERPMG_GETWEIGHTS(
     & nref
     & ,ibreflo0,ibreflo1,ibreflo2
     & ,ibrefhi0,ibrefhi1,ibrefhi2
     & ,wtcrnr
     & ,iwtcrnrlo0,iwtcrnrlo1,iwtcrnrlo2
     & ,iwtcrnrhi0,iwtcrnrhi1,iwtcrnrhi2
     & ,nwtcrnrcomp
     & )
      implicit none
      integer nref
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer nwtcrnrcomp
      integer iwtcrnrlo0,iwtcrnrlo1,iwtcrnrlo2
      integer iwtcrnrhi0,iwtcrnrhi1,iwtcrnrhi2
      REAL*8 wtcrnr(
     & iwtcrnrlo0:iwtcrnrhi0,
     & iwtcrnrlo1:iwtcrnrhi1,
     & iwtcrnrlo2:iwtcrnrhi2,
     & 0:nwtcrnrcomp-1)
      integer iref0,iref1,iref2
      integer ib0,ib1,ib2
      integer nvwt
      integer ibmax0,ibmax1,ibmax2
      REAL*8 refinv, wt
      REAL*8 fraci0,fraci1,fraci2
      REAL*8 wti0,wti1,wti2
      refinv = (1.0d0) / nref
      nvwt = 0
      do iref2 = ibreflo2,ibrefhi2
      do iref1 = ibreflo1,ibrefhi1
      do iref0 = ibreflo0,ibrefhi0
         call maxb(iref0, ibmax0)
         call maxb(iref1, ibmax1)
         call maxb(iref2, ibmax2)
         fraci0 = iref0 * refinv
         fraci1 = iref1 * refinv
         fraci2 = iref2 * refinv
         do ib0 = 0, ibmax0
         do ib1 = 0, ibmax1
         do ib2 = 0, ibmax2
            wt = (1.0d0)
            call wtside(ib0, fraci0, wti0)
            wt = wt * wti0
            call wtside(ib1, fraci1, wti1)
            wt = wt * wti1
            call wtside(ib2, fraci2, wti2)
            wt = wt * wti2
            wtcrnr( ib0,ib1,ib2, nvwt ) = wt
         end do
         end do
         end do
         nvwt = nvwt + 1
      enddo
      enddo
      enddo
      return
      end
      subroutine NODEINTERPMG(
     & fine
     & ,ifinelo0,ifinelo1,ifinelo2
     & ,ifinehi0,ifinehi1,ifinehi2
     & ,nfinecomp
     & ,crse
     & ,icrselo0,icrselo1,icrselo2
     & ,icrsehi0,icrsehi1,icrsehi2
     & ,ncrsecomp
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,nref
     & ,ibreflo0,ibreflo1,ibreflo2
     & ,ibrefhi0,ibrefhi1,ibrefhi2
     & ,wtcrnr
     & ,iwtcrnrlo0,iwtcrnrlo1,iwtcrnrlo2
     & ,iwtcrnrhi0,iwtcrnrhi1,iwtcrnrhi2
     & ,nwtcrnrcomp
     & )
      implicit none
      integer nfinecomp
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL*8 fine(
     & ifinelo0:ifinehi0,
     & ifinelo1:ifinehi1,
     & ifinelo2:ifinehi2,
     & 0:nfinecomp-1)
      integer ncrsecomp
      integer icrselo0,icrselo1,icrselo2
      integer icrsehi0,icrsehi1,icrsehi2
      REAL*8 crse(
     & icrselo0:icrsehi0,
     & icrselo1:icrsehi1,
     & icrselo2:icrsehi2,
     & 0:ncrsecomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer nref
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer nwtcrnrcomp
      integer iwtcrnrlo0,iwtcrnrlo1,iwtcrnrlo2
      integer iwtcrnrhi0,iwtcrnrhi1,iwtcrnrhi2
      REAL*8 wtcrnr(
     & iwtcrnrlo0:iwtcrnrhi0,
     & iwtcrnrlo1:iwtcrnrhi1,
     & iwtcrnrlo2:iwtcrnrhi2,
     & 0:nwtcrnrcomp-1)
      integer iref0,iref1,iref2, icrse0,icrse1,icrse2
      integer ifine0,ifine1,ifine2, ib0,ib1,ib2;
      integer var, ncomp, nvwt
      integer ibmax0,ibmax1,ibmax2
      integer icmin0,icmin1,icmin2
      integer icmax0,icmax1,icmax2;
      REAL*8 csum, refinv
      ncomp = nfinecomp
      if (ncomp .ne. ncrsecomp) then
         print *, 'nodeinterpmg: fine and crse incompatible'
         call MAYDAY_ERROR()
      endif
      refinv = (1.0d0) / nref
      icmin0 = iregionlo0
      icmin1 = iregionlo1
      icmin2 = iregionlo2
      nvwt = 0
      do iref2 = ibreflo2,ibrefhi2
      do iref1 = ibreflo1,ibrefhi1
      do iref0 = ibreflo0,ibrefhi0
         call maxb(iref0, ibmax0)
         call maxb(iref1, ibmax1)
         call maxb(iref2, ibmax2)
         icmax0 = iregionhi0 + (1-ibmax0)
         icmax1 = iregionhi1 + (1-ibmax1)
         icmax2 = iregionhi2 + (1-ibmax2)
         do icrse0 = icmin0, icmax0
         do icrse1 = icmin1, icmax1
         do icrse2 = icmin2, icmax2
            ifine0 = nref*icrse0 + iref0
            ifine1 = nref*icrse1 + iref1
            ifine2 = nref*icrse2 + iref2
            do var = 0, ncomp-1
               csum = 0
               do ib0 = 0, ibmax0
               do ib1 = 0, ibmax1
               do ib2 = 0, ibmax2
                  csum = csum + wtcrnr( ib0,ib1,ib2, nvwt ) *
     & crse( icrse0+ib0,icrse1+ib1,icrse2+ib2, var)
               end do
               end do
               end do
               fine ( ifine0,ifine1,ifine2, var ) = csum +
     & fine ( ifine0,ifine1,ifine2, var )
            end do
         end do
         end do
         end do
         nvwt = nvwt + 1
      enddo
      enddo
      enddo
      return
      end
      subroutine WTSIDE(
     & i
     & ,frac
     & ,wt
     & )
      implicit none
      integer i
      REAL*8 frac
      REAL*8 wt
      if (i .eq. 0) then
         wt = (1.0d0) - frac
      else
         wt = frac
      endif
      return
      end
      subroutine MAXB(
     & iref
     & ,ibmax
     & )
      implicit none
      integer iref
      integer ibmax
      if (iref .eq. 0) then
         ibmax = 0
      else
         ibmax = 1
      endif
      return
      end
