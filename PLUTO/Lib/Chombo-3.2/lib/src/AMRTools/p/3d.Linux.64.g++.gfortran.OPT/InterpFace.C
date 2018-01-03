#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine INTERPFACECONSTANT(
     &           fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,nfinecomp
     &           ,coarse
     &           ,icoarselo0,icoarselo1,icoarselo2
     &           ,icoarsehi0,icoarsehi1,icoarsehi2
     &           ,ncoarsecomp
     &           ,iblo0,iblo1,iblo2
     &           ,ibhi0,ibhi1,ibhi2
     &           ,ref_ratio
     &           ,ibreflo0,ibreflo1,ibreflo2
     &           ,ibrefhi0,ibrefhi1,ibrefhi2
     &           ,dir
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
      integer ncoarsecomp
      integer icoarselo0,icoarselo1,icoarselo2
      integer icoarsehi0,icoarsehi1,icoarsehi2
      REAL_T coarse(
     &           icoarselo0:icoarsehi0,
     &           icoarselo1:icoarsehi1,
     &           icoarselo2:icoarsehi2,
     &           0:ncoarsecomp-1)
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer ref_ratio
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer dir
      integer var
      integer ic0,ic1,ic2
      integer ifine0,ifine1,ifine2
      integer ii0,ii1,ii2
      do var = 0, ncoarsecomp - 1
         
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

            
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

            
               ifine0 = ic0*ref_ratio + ii0
               ifine1 = ic1*ref_ratio + ii1
               ifine2 = ic2*ref_ratio + ii2
               fine(ifine0,ifine1,ifine2,var) = coarse(ic0,ic1,ic2,var)
            
      enddo
      enddo
      enddo
         
      enddo
      enddo
      enddo
      end do
      return
      end
      subroutine INTERPLINEARFACE(
     &           fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,nfinecomp
     &           ,slope
     &           ,islopelo0,islopelo1,islopelo2
     &           ,islopehi0,islopehi1,islopehi2
     &           ,nslopecomp
     &           ,iblo0,iblo1,iblo2
     &           ,ibhi0,ibhi1,ibhi2
     &           ,dir
     &           ,ref_ratio
     &           ,ibreffacelo0,ibreffacelo1,ibreffacelo2
     &           ,ibreffacehi0,ibreffacehi1,ibreffacehi2
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
      integer nslopecomp
      integer islopelo0,islopelo1,islopelo2
      integer islopehi0,islopehi1,islopehi2
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1,
     &           islopelo2:islopehi2,
     &           0:nslopecomp-1)
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer dir
      integer ref_ratio
      integer ibreffacelo0,ibreffacelo1,ibreffacelo2
      integer ibreffacehi0,ibreffacehi1,ibreffacehi2
      integer  ic 0, ic 1, ic 2
      integer  ifine 0, ifine 1, ifine 2
      integer  ii 0, ii 1, ii 2
      integer var, id
      REAL_T dxf
      do var = 0, nfinecomp - 1
         
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

              
      do ii2 = ibreffacelo2,ibreffacehi2
      do ii1 = ibreffacelo1,ibreffacehi1
      do ii0 = ibreffacelo0,ibreffacehi0

              
                  ifine0 = ic0*ref_ratio + ii0
                  ifine1 = ic1*ref_ratio + ii1
                  ifine2 = ic2*ref_ratio + ii2
              
                  if (dir .eq. 0) then
                      id = ii0
                  else if (dir .eq. 1) then
                      id = ii1
                  else if (dir .eq. 2) then
                      id = ii2
                  endif
                  dxf = -half + ( (id+half) / ref_ratio )
                  fine( ifine0,ifine1,ifine2,var) =
     &                 fine( ifine0,ifine1,ifine2,var) +
     &                 dxf * slope (   ic 0,  ic 1,  ic 2, var )
              
      enddo
      enddo
      enddo
          
      enddo
      enddo
      enddo
      end do
      return
      end
      subroutine INTERPLINEARINTERIORFACE(
     &           fine
     &           ,ifinelo0,ifinelo1,ifinelo2
     &           ,ifinehi0,ifinehi1,ifinehi2
     &           ,nfinecomp
     &           ,ibcoarselo0,ibcoarselo1,ibcoarselo2
     &           ,ibcoarsehi0,ibcoarsehi1,ibcoarsehi2
     &           ,ref_ratio
     &           ,facedir
     &           ,iinteriorrefboxlo0,iinteriorrefboxlo1,iinteriorrefboxlo2
     &           ,iinteriorrefboxhi0,iinteriorrefboxhi1,iinteriorrefboxhi2
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nfinecomp
      integer ifinelo0,ifinelo1,ifinelo2
      integer ifinehi0,ifinehi1,ifinehi2
      REAL_T fine(
     &           ifinelo0:ifinehi0,
     &           ifinelo1:ifinehi1,
     &           ifinelo2:ifinehi2,
     &           0:nfinecomp-1)
      integer ibcoarselo0,ibcoarselo1,ibcoarselo2
      integer ibcoarsehi0,ibcoarsehi1,ibcoarsehi2
      integer ref_ratio
      integer facedir
      integer iinteriorrefboxlo0,iinteriorrefboxlo1,iinteriorrefboxlo2
      integer iinteriorrefboxhi0,iinteriorrefboxhi1,iinteriorrefboxhi2
      integer ic0,ic1,ic2
      integer ifine0,ifine1,ifine2
      integer ii0,ii1,ii2
      integer iloface0,iloface1,iloface2
      integer ihiface0,ihiface1,ihiface2
      integer var, id
      REAL_T dxf, diff
      REAL_T loval, hival
      do var=0, nfinecomp -1
         
      do ic2 = ibcoarselo2,ibcoarsehi2
      do ic1 = ibcoarselo1,ibcoarsehi1
      do ic0 = ibcoarselo0,ibcoarsehi0

            
      do ii2 = iinteriorrefboxlo2,iinteriorrefboxhi2
      do ii1 = iinteriorrefboxlo1,iinteriorrefboxhi1
      do ii0 = iinteriorrefboxlo0,iinteriorrefboxhi0

            
              ifine0 = ic0*ref_ratio + ii0
              ifine1 = ic1*ref_ratio + ii1
              ifine2 = ic2*ref_ratio + ii2
              
              iloface0 = ic0*ref_ratio + (1-CHF_ID(0,facedir))*ii0
              iloface1 = ic1*ref_ratio + (1-CHF_ID(1,facedir))*ii1
              iloface2 = ic2*ref_ratio + (1-CHF_ID(2,facedir))*ii2
              
              ihiface0 = iloface0 + ref_ratio*CHF_ID(0,facedir)
              ihiface1 = iloface1 + ref_ratio*CHF_ID(1,facedir)
              ihiface2 = iloface2 + ref_ratio*CHF_ID(2,facedir)
              
              if (facedir .eq. 0) then
                 id = ii0
              else if (facedir .eq. 1) then
                 id = ii1
              else if (facedir .eq. 2) then
                 id = ii2
              endif
              dxf = float(id)/ref_ratio
              diff = fine(ihiface0,ihiface1,ihiface2,var)
     &                -fine(iloface0,iloface1,iloface2,var)
              fine( ifine0,ifine1,ifine2,var) =
     &            fine(iloface0,iloface1,iloface2,var)
     &           +dxf*diff
            
      enddo
      enddo
      enddo
          
      enddo
      enddo
      enddo
       enddo
       return
       end
      subroutine INTERPLIMITFACE(
     &           islope
     &           ,iislopelo0,iislopelo1,iislopelo2
     &           ,iislopehi0,iislopehi1,iislopehi2
     &           ,nislopecomp
     &           ,jslope
     &           ,ijslopelo0,ijslopelo1,ijslopelo2
     &           ,ijslopehi0,ijslopehi1,ijslopehi2
     &           ,njslopecomp
     &           ,kslope
     &           ,ikslopelo0,ikslopelo1,ikslopelo2
     &           ,ikslopehi0,ikslopehi1,ikslopehi2
     &           ,nkslopecomp
     &           ,lslope
     &           ,ilslopelo0,ilslopelo1,ilslopelo2
     &           ,ilslopehi0,ilslopehi1,ilslopehi2
     &           ,nlslopecomp
     &           ,mslope
     &           ,imslopelo0,imslopelo1,imslopelo2
     &           ,imslopehi0,imslopehi1,imslopehi2
     &           ,nmslopecomp
     &           ,nslope
     &           ,inslopelo0,inslopelo1,inslopelo2
     &           ,inslopehi0,inslopehi1,inslopehi2
     &           ,nnslopecomp
     &           ,state
     &           ,istatelo0,istatelo1,istatelo2
     &           ,istatehi0,istatehi1,istatehi2
     &           ,nstatecomp
     &           ,iblo0,iblo1,iblo2
     &           ,ibhi0,ibhi1,ibhi2
     &           ,ibnlo0,ibnlo1,ibnlo2
     &           ,ibnhi0,ibnhi1,ibnhi2
     &           ,ivalidBoxlo0,ivalidBoxlo1,ivalidBoxlo2
     &           ,ivalidBoxhi0,ivalidBoxhi1,ivalidBoxhi2
     &           ,normaldir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nislopecomp
      integer iislopelo0,iislopelo1,iislopelo2
      integer iislopehi0,iislopehi1,iislopehi2
      REAL_T islope(
     &           iislopelo0:iislopehi0,
     &           iislopelo1:iislopehi1,
     &           iislopelo2:iislopehi2,
     &           0:nislopecomp-1)
      integer njslopecomp
      integer ijslopelo0,ijslopelo1,ijslopelo2
      integer ijslopehi0,ijslopehi1,ijslopehi2
      REAL_T jslope(
     &           ijslopelo0:ijslopehi0,
     &           ijslopelo1:ijslopehi1,
     &           ijslopelo2:ijslopehi2,
     &           0:njslopecomp-1)
      integer nkslopecomp
      integer ikslopelo0,ikslopelo1,ikslopelo2
      integer ikslopehi0,ikslopehi1,ikslopehi2
      REAL_T kslope(
     &           ikslopelo0:ikslopehi0,
     &           ikslopelo1:ikslopehi1,
     &           ikslopelo2:ikslopehi2,
     &           0:nkslopecomp-1)
      integer nlslopecomp
      integer ilslopelo0,ilslopelo1,ilslopelo2
      integer ilslopehi0,ilslopehi1,ilslopehi2
      REAL_T lslope(
     &           ilslopelo0:ilslopehi0,
     &           ilslopelo1:ilslopehi1,
     &           ilslopelo2:ilslopehi2,
     &           0:nlslopecomp-1)
      integer nmslopecomp
      integer imslopelo0,imslopelo1,imslopelo2
      integer imslopehi0,imslopehi1,imslopehi2
      REAL_T mslope(
     &           imslopelo0:imslopehi0,
     &           imslopelo1:imslopehi1,
     &           imslopelo2:imslopehi2,
     &           0:nmslopecomp-1)
      integer nnslopecomp
      integer inslopelo0,inslopelo1,inslopelo2
      integer inslopehi0,inslopehi1,inslopehi2
      REAL_T nslope(
     &           inslopelo0:inslopehi0,
     &           inslopelo1:inslopehi1,
     &           inslopelo2:inslopehi2,
     &           0:nnslopecomp-1)
      integer nstatecomp
      integer istatelo0,istatelo1,istatelo2
      integer istatehi0,istatehi1,istatehi2
      REAL_T state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1,
     &           istatelo2:istatehi2,
     &           0:nstatecomp-1)
      integer iblo0,iblo1,iblo2
      integer ibhi0,ibhi1,ibhi2
      integer ibnlo0,ibnlo1,ibnlo2
      integer ibnhi0,ibnhi1,ibnhi2
      integer ivalidBoxlo0,ivalidBoxlo1,ivalidBoxlo2
      integer ivalidBoxhi0,ivalidBoxhi1,ivalidBoxhi2
      integer normaldir
      integer   i 0,  i 1,  i 2, var
      integer   ii 0,  ii 1,  ii 2
      integer   in 0,  in 1,  in 2
      REAL_T statemax, statemin, deltasum, etamax, etamin, eta
      REAL_T tempone, tempzero
      tempone = one
      tempzero = 0
      do var = 0, nislopecomp - 1
         
      do i2 = iblo2,ibhi2
      do i1 = iblo1,ibhi1
      do i0 = iblo0,ibhi0

             statemax = state ( i0,i1,i2, var )
             statemin = state ( i0,i1,i2, var )
             
      do ii2 = ibnlo2,ibnhi2
      do ii1 = ibnlo1,ibnhi1
      do ii0 = ibnlo0,ibnhi0

             
                 in0 = i0 + ii0
                 in1 = i1 + ii1
                 in2 = i2 + ii2
                 if (
                 
     &                in0 .ge. ivalidBoxlo0 .and.
     &                in0 .le. ivalidBoxhi0 
     &                .and.
     &                in1 .ge. ivalidBoxlo1 .and.
     &                in1 .le. ivalidBoxhi1 
     &                .and.
     &                in2 .ge. ivalidBoxlo2 .and.
     &                in2 .le. ivalidBoxhi2 
     &                )
     &        then
                    statemax = max ( statemax, state(in0,in1,in2,var))
                    statemin = min ( statemin, state(in0,in1,in2,var))
                 endif
             
      enddo
      enddo
      enddo
             deltasum = half * (
                
     &            (1-CHF_ID(normaldir,0))*abs(islope(i0,i1,i2,var))
     &            +
     &            (1-CHF_ID(normaldir,1))*abs(jslope(i0,i1,i2,var))
     &            +
     &            (1-CHF_ID(normaldir,2))*abs(kslope(i0,i1,i2,var))
     &            )
             if ( deltasum .gt. 0 ) then
                etamax = ( statemax - state ( i0,i1,i2, var ) )
     &               / deltasum
                etamin = ( state ( i0,i1,i2, var ) - statemin )
     &               / deltasum
                eta = max ( min ( etamin, etamax, tempone ), tempzero )
                
                islope ( i0,i1,i2, var ) =
     &               eta * islope ( i0,i1,i2, var ) 
                jslope ( i0,i1,i2, var ) =
     &               eta * jslope ( i0,i1,i2, var ) 
                kslope ( i0,i1,i2, var ) =
     &               eta * kslope ( i0,i1,i2, var ) 
             end if
         
      enddo
      enddo
      enddo
      end do
      return
      end
