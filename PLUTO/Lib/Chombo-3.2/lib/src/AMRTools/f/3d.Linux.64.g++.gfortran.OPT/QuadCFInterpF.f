      subroutine QUADINTERP(
     & phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,nphicomp
     & ,phistar
     & ,iphistarlo0,iphistarlo1,iphistarlo2
     & ,iphistarhi0,iphistarhi1,iphistarhi2
     & ,nphistarcomp
     & ,iboxlo0,iboxlo1,iboxlo2
     & ,iboxhi0,iboxhi1,iboxhi2
     & ,ihilo
     & ,h
     & ,idir
     & ,scomp
     & ,ecomp
     & ,nref
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2,
     & 0:nphicomp-1)
      integer nphistarcomp
      integer iphistarlo0,iphistarlo1,iphistarlo2
      integer iphistarhi0,iphistarhi1,iphistarhi2
      REAL*8 phistar(
     & iphistarlo0:iphistarhi0,
     & iphistarlo1:iphistarhi1,
     & iphistarlo2:iphistarhi2,
     & 0:nphistarcomp-1)
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer ihilo
      REAL*8 h
      integer idir
      integer scomp
      integer ecomp
      integer nref
      integer i0,i1,i2
      integer ii0,ii1,ii2
      integer n
      REAL*8 x, pa, pb, ps, a, b, frac, denom, xsquared
      REAL*8 mult, invh
      frac = (2.0d0)/(h*h)
      denom = nref*nref + 4*nref + 3
      mult = frac/denom
      invh = (1.0d0) / h
      x = (2.0d0) * h
      xsquared = (4.0d0) * h*h
      ii0= ihilo*CHF_ID(0, idir)
      ii1= ihilo*CHF_ID(1, idir)
      ii2= ihilo*CHF_ID(2, idir)
      do n=scomp,ecomp
      do i2 = iboxlo2,iboxhi2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0
            pa = phi(i0 -2*ii0,i1 -2*ii1,i2 -2*ii2,n)
            pb = phi(i0 -ii0,i1 -ii1,i2 -ii2,n)
            ps = phistar(i0 +ii0,i1 +ii1,i2 +ii2,n)
            a = mult*((2.0d0)*ps + (nref+1)*pa - (nref+3)*pb)
            b = (pb-pa)*invh - a*h
            phi(i0,i1,i2,n) = xsquared*a + b*x + pa
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine PHISTAR(
     & fPhiStar
     & ,ifPhiStarlo0,ifPhiStarlo1,ifPhiStarlo2
     & ,ifPhiStarhi0,ifPhiStarhi1,ifPhiStarhi2
     & ,nfPhiStarcomp
     & ,iregionlo0,iregionlo1,iregionlo2
     & ,iregionhi0,iregionhi1,iregionhi2
     & ,phic
     & ,iphiclo0,iphiclo1,iphiclo2
     & ,iphichi0,iphichi1,iphichi2
     & ,nphiccomp
     & ,coarslope
     & ,icoarslopelo0,icoarslopelo1,icoarslopelo2
     & ,icoarslopehi0,icoarslopehi1,icoarslopehi2
     & ,ncoarslopecomp
     & ,coarcurva
     & ,icoarcurvalo0,icoarcurvalo1,icoarcurvalo2
     & ,icoarcurvahi0,icoarcurvahi1,icoarcurvahi2
     & ,ncoarcurvacomp
     & ,coarmixed
     & ,icoarmixedlo0,icoarmixedlo1,icoarmixedlo2
     & ,icoarmixedhi0,icoarmixedhi1,icoarmixedhi2
     & ,ncoarmixedcomp
     & ,dxf
     & ,ivar
     & ,dir
     & ,sign
     & ,nRef
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer nfPhiStarcomp
      integer ifPhiStarlo0,ifPhiStarlo1,ifPhiStarlo2
      integer ifPhiStarhi0,ifPhiStarhi1,ifPhiStarhi2
      REAL*8 fPhiStar(
     & ifPhiStarlo0:ifPhiStarhi0,
     & ifPhiStarlo1:ifPhiStarhi1,
     & ifPhiStarlo2:ifPhiStarhi2,
     & 0:nfPhiStarcomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      integer nphiccomp
      integer iphiclo0,iphiclo1,iphiclo2
      integer iphichi0,iphichi1,iphichi2
      REAL*8 phic(
     & iphiclo0:iphichi0,
     & iphiclo1:iphichi1,
     & iphiclo2:iphichi2,
     & 0:nphiccomp-1)
      integer ncoarslopecomp
      integer icoarslopelo0,icoarslopelo1,icoarslopelo2
      integer icoarslopehi0,icoarslopehi1,icoarslopehi2
      REAL*8 coarslope(
     & icoarslopelo0:icoarslopehi0,
     & icoarslopelo1:icoarslopehi1,
     & icoarslopelo2:icoarslopehi2,
     & 0:ncoarslopecomp-1)
      integer ncoarcurvacomp
      integer icoarcurvalo0,icoarcurvalo1,icoarcurvalo2
      integer icoarcurvahi0,icoarcurvahi1,icoarcurvahi2
      REAL*8 coarcurva(
     & icoarcurvalo0:icoarcurvahi0,
     & icoarcurvalo1:icoarcurvahi1,
     & icoarcurvalo2:icoarcurvahi2,
     & 0:ncoarcurvacomp-1)
      integer ncoarmixedcomp
      integer icoarmixedlo0,icoarmixedlo1,icoarmixedlo2
      integer icoarmixedhi0,icoarmixedhi1,icoarmixedhi2
      REAL*8 coarmixed(
     & icoarmixedlo0:icoarmixedhi0,
     & icoarmixedlo1:icoarmixedhi1,
     & icoarmixedlo2:icoarmixedhi2,
     & 0:ncoarmixedcomp-1)
      REAL*8 dxf
      integer ivar
      integer dir
      integer sign
      integer nRef
      REAL*8 xf1, xc1, xf2, xc2, x1, x2, dxc
      REAL*8 aa, update1, update2, update3
      integer i0,i1,i2
      integer ii0,ii1,ii2
      integer ir0,ir1,ir2
      integer ic(0:3 -1)
      integer ivf(0:3 -1)
      integer YOU(1:2, 0:2), you1, you2
      data YOU / 1, 2, 0, 2, 0, 1 /
      dxc = nRef * dxf
      you1 = YOU(1,dir)
      you2 = YOU(2,dir)
      ii0= sign*CHF_ID(0, dir)
      ii1= sign*CHF_ID(1, dir)
      ii2= sign*CHF_ID(2, dir)
      do ir2 = iregionlo2,iregionhi2
      do ir1 = iregionlo1,iregionhi1
      do ir0 = iregionlo0,iregionhi0
         ic(0)=ir0/nRef
         ic(1)=ir1/nRef
         ic(2)=ir2/nRef
         ivf(0)=ir0
         ivf(1)=ir1
         ivf(2)=ir2
         i0=ir0+ii0
         i1=ir1+ii1
         i2=ir2+ii2
         xf1 = (ivf(you1)+(0.500d0))*dxf
         xc1 = ( ic(you1)+(0.500d0))*dxc
         xf2 = (ivf(you2)+(0.500d0))*dxf
         xc2 = ( ic(you2)+(0.500d0))*dxc
         x1 = xf1-xc1
         x2 = xf2-xc2
         aa= phic(ic(0),ic(1),ic(2),ivar)
         update1=x1*coarslope(ic(0),ic(1),ic(2),you1) +
     & (0.500d0)*x1*x1*coarcurva(ic(0),ic(1),ic(2),you1)
         update2=x2*coarslope(ic(0),ic(1),ic(2),you2) +
     & (0.500d0)*x2*x2*coarcurva(ic(0),ic(1),ic(2),you2)
         update3=x1*x2*coarmixed(ic(0),ic(1),ic(2),0)
         fPhiStar(i0,i1,i2,ivar) = aa+update1+update2+update3
      enddo
      enddo
      enddo
      return
      end
