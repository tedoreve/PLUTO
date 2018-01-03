#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine FACENODEBC(
     &           state
     &           ,istatelo0,istatelo1,istatelo2
     &           ,istatehi0,istatehi1,istatehi2
     &           ,nstatecomp
     &           ,neumfac
     &           ,ineumfaclo0,ineumfaclo1,ineumfaclo2
     &           ,ineumfachi0,ineumfachi1,ineumfachi2
     &           ,nneumfaccomp
     &           ,dircfac
     &           ,idircfaclo0,idircfaclo1,idircfaclo2
     &           ,idircfachi0,idircfachi1,idircfachi2
     &           ,ndircfaccomp
     &           ,inhmval
     &           ,iinhmvallo0,iinhmvallo1,iinhmvallo2
     &           ,iinhmvalhi0,iinhmvalhi1,iinhmvalhi2
     &           ,ninhmvalcomp
     &           ,ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
     &           ,ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
     &           ,idir
     &           ,side
     &           ,dx
     &           ,startcomp
     &           ,endcomp
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nstatecomp
      integer istatelo0,istatelo1,istatelo2
      integer istatehi0,istatehi1,istatehi2
      REAL_T state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1,
     &           istatelo2:istatehi2,
     &           0:nstatecomp-1)
      integer nneumfaccomp
      integer ineumfaclo0,ineumfaclo1,ineumfaclo2
      integer ineumfachi0,ineumfachi1,ineumfachi2
      REAL_T neumfac(
     &           ineumfaclo0:ineumfachi0,
     &           ineumfaclo1:ineumfachi1,
     &           ineumfaclo2:ineumfachi2,
     &           0:nneumfaccomp-1)
      integer ndircfaccomp
      integer idircfaclo0,idircfaclo1,idircfaclo2
      integer idircfachi0,idircfachi1,idircfachi2
      REAL_T dircfac(
     &           idircfaclo0:idircfachi0,
     &           idircfaclo1:idircfachi1,
     &           idircfaclo2:idircfachi2,
     &           0:ndircfaccomp-1)
      integer ninhmvalcomp
      integer iinhmvallo0,iinhmvallo1,iinhmvallo2
      integer iinhmvalhi0,iinhmvalhi1,iinhmvalhi2
      REAL_T inhmval(
     &           iinhmvallo0:iinhmvalhi0,
     &           iinhmvallo1:iinhmvalhi1,
     &           iinhmvallo2:iinhmvalhi2,
     &           0:ninhmvalcomp-1)
      integer ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
      integer ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
      integer idir
      integer side
      REAL_T dx
      integer startcomp
      integer endcomp
      REAL_T nfac, dfac, ival, sval,denom,numer
      integer ncomp,nc
      integer i0,i1,i2, ii0,ii1,ii2
      ncomp = nstatecomp
      if(ncomp .ne. nneumfaccomp) then
          call MAYDAY_ERROR()
      endif
      if(ncomp .ne. ndircfaccomp) then
          call MAYDAY_ERROR()
      endif
      if(ncomp .ne. ninhmvalcomp) then
          call MAYDAY_ERROR()
      endif
      if ((side .ne. -1) .and. (side .ne. 1)) then
          call MAYDAY_ERROR()
      endif
      
      ii0= side*CHF_ID(0, idir)

      ii1= side*CHF_ID(1, idir)

      ii2= side*CHF_ID(2, idir)

      do nc = startcomp, endcomp
          
      do i2 = ifaceboxlo2,ifaceboxhi2
      do i1 = ifaceboxlo1,ifaceboxhi1
      do i0 = ifaceboxlo0,ifaceboxhi0

              nfac = neumfac(i0,i1,i2, nc)
              dfac = dircfac(i0,i1,i2, nc)
              ival = inhmval(i0,i1,i2, nc)
              sval = state(i0-ii0,i1-ii1,i2-ii2, nc)
              denom = dfac + side*(nfac/dx)
              numer = ival + side*(nfac/dx)*sval
#ifndef NDEBUG
              if (abs(denom) .lt. 1.0e-9) then
                  call MAYDAY_ERROR()
              endif
#endif
              state(i0,i1,i2, nc) = numer/denom
          
      enddo
      enddo
      enddo
      enddo
      return
      end
