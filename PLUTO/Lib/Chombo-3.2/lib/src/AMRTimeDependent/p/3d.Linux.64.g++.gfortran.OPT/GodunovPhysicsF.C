#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine FLUXDIFFF(
     &           diff
     &           ,idifflo0,idifflo1,idifflo2
     &           ,idiffhi0,idiffhi1,idiffhi2
     &           ,ndiffcomp
     &           ,F
     &           ,iFlo0,iFlo1,iFlo2
     &           ,iFhi0,iFhi1,iFhi2
     &           ,nFcomp
     &           ,idir
     &           ,iboxlo0,iboxlo1,iboxlo2
     &           ,iboxhi0,iboxhi1,iboxhi2
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer ndiffcomp
      integer idifflo0,idifflo1,idifflo2
      integer idiffhi0,idiffhi1,idiffhi2
      REAL_T diff(
     &           idifflo0:idiffhi0,
     &           idifflo1:idiffhi1,
     &           idifflo2:idiffhi2,
     &           0:ndiffcomp-1)
      integer nFcomp
      integer iFlo0,iFlo1,iFlo2
      integer iFhi0,iFhi1,iFhi2
      REAL_T F(
     &           iFlo0:iFhi0,
     &           iFlo1:iFhi1,
     &           iFlo2:iFhi2,
     &           0:nFcomp-1)
      integer idir
      integer iboxlo0,iboxlo1,iboxlo2
      integer iboxhi0,iboxhi1,iboxhi2
      integer i0,i1,i2
      integer c2fLo0,c2fLo1,c2fLo2
      integer c2fHi0,c2fHi1,c2fHi2
      integer iv
      
      c2fLo0= 0*CHF_ID(0, idir)

      c2fLo1= 0*CHF_ID(1, idir)

      c2fLo2= 0*CHF_ID(2, idir)

      
      c2fHi0= 1*CHF_ID(0, idir)

      c2fHi1= 1*CHF_ID(1, idir)

      c2fHi2= 1*CHF_ID(2, idir)

      do iv = 0,ndiffcomp - 1
         
      do i2 = iboxlo2,iboxhi2
      do i1 = iboxlo1,iboxhi1
      do i0 = iboxlo0,iboxhi0

            diff(i0,i1,i2, iv) =
     &        F(i0 +c2fHi0,i1 +c2fHi1,i2 +c2fHi2, iv) -
     &        F(i0 +c2fLo0,i1 +c2fLo1,i2 +c2fLo2, iv)
         
      enddo
      enddo
      enddo
      enddo
      return
      end
