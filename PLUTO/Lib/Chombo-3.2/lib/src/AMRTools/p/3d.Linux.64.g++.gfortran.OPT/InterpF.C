#include "REAL.H"
#include "SPACE.H"
#include "CONSTANTS.H"

#include "CONSTANTS.H"
      subroutine INTERPCONSTANT(
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
     &           ,dx
     &           ,stretch
     &           ,xbeg
     &           ,ixbeghi0
     &           ,ibreflo0,ibreflo1,ibreflo2
     &           ,ibrefhi0,ibrefhi1,ibrefhi2
     &           ,geometry
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
      REAL_T dx
      REAL_T stretch
      integer ixbeghi0
      REAL_T xbeg(
     &           0:ixbeghi0)
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer geometry
      integer var
      integer ic0,ic1,ic2
      integer if0,if1,if2
      integer ii0,ii1,ii2
      REAL_T  xlf0,xlf1,xlf2
      REAL_T  xrf0,xrf1,xrf2
      REAL_T  volume
      if (geometry .eq. 1) then
      do var = 0, ncoarsecomp - 1
         
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

            
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

            
               if0 = ic0*ref_ratio + ii0
               if1 = ic1*ref_ratio + ii1
               if2 = ic2*ref_ratio + ii2
               fine(if0,if1,if2,var) = coarse(ic0,ic1,ic2,var)
            
      enddo
      enddo
      enddo
         
      enddo
      enddo
      enddo
      end do
      endif
      if ((geometry .eq. 2) .or. (geometry .eq. 5)) then
      do var = 0, ncoarsecomp - 1
         
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

            
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

            
               if0 = ic0*ref_ratio + ii0
               volume = abs(xbeg(0)/dx+if0+half)
               if1 = ic1*ref_ratio + ii1
               if2 = ic2*ref_ratio + ii2
               fine(if0,if1,if2,var) = coarse(ic0,ic1,ic2,var)*volume
            
      enddo
      enddo
      enddo
         
      enddo
      enddo
      enddo
      end do
      endif
      if (geometry .eq. 3) then
      do var = 0, ncoarsecomp - 1
         
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

            
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

            
               if0 = ic0*ref_ratio + ii0
               xlf0 = xbeg(0)+if0*dx
               xrf0 = xlf0 + dx
               volume = (xrf0*xrf0*xrf0-xlf0*xlf0*xlf0)*third
               if1 = ic1*ref_ratio + ii1
               xlf1 = xbeg(1)+if1*dx*stretch
               xrf1 = xlf1 + dx*stretch
               volume = volume*(cos(xlf1)-cos(xrf1))
               if2 = ic2*ref_ratio + ii2
               fine(if0,if1,if2,var) = coarse(ic0,ic1,ic2,var)*volume
            
      enddo
      enddo
      enddo
         
      enddo
      enddo
      enddo
      end do
      endif
      if (geometry .eq. 4) then
      do var = 0, ncoarsecomp - 1
         
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

            
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

            
               if0 = ic0*ref_ratio + ii0
               xlf0 = if0*dx
               xrf0 = xlf0 + dx
               volume = (exp(three*xrf0)-exp(three*xlf0))*third
               if1 = ic1*ref_ratio + ii1
               xlf1 = xbeg(1)+if1*dx*stretch
               xrf1 = xlf1 + dx*stretch
               volume = volume*(cos(xlf1)-cos(xrf1))
               if2 = ic2*ref_ratio + ii2
               fine(if0,if1,if2,var) = coarse(ic0,ic1,ic2,var)*volume
            
      enddo
      enddo
      enddo
         
      enddo
      enddo
      enddo
      end do
      endif
      if (geometry .eq. 6) then
      do var = 0, ncoarsecomp - 1
         
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

            
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

            
               if0 = ic0*ref_ratio + ii0
               xlf0 = if0*dx
               xrf0 = xlf0 + dx
               volume = (exp(two*xrf0)-exp(two*xlf0))*half
               if1 = ic1*ref_ratio + ii1
               if2 = ic2*ref_ratio + ii2
               fine(if0,if1,if2,var) = coarse(ic0,ic1,ic2,var)*volume
            
      enddo
      enddo
      enddo
         
      enddo
      enddo
      enddo
      end do
      endif
      return
      end
      subroutine INTERPCENTRALSLOPE(
     &           slope
     &           ,islopelo0,islopelo1,islopelo2
     &           ,islopehi0,islopehi1,islopehi2
     &           ,nslopecomp
     &           ,state
     &           ,istatelo0,istatelo1,istatelo2
     &           ,istatehi0,istatehi1,istatehi2
     &           ,nstatecomp
     &           ,iblo0,iblo1,iblo2
     &           ,ibhi0,ibhi1,ibhi2
     &           ,dir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nslopecomp
      integer islopelo0,islopelo1,islopelo2
      integer islopehi0,islopehi1,islopehi2
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1,
     &           islopelo2:islopehi2,
     &           0:nslopecomp-1)
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
      integer dir
      integer i0,i1,i2
      integer ii0,ii1,ii2
      integer var
      
      ii0=CHF_ID(0, dir)

      ii1=CHF_ID(1, dir)

      ii2=CHF_ID(2, dir)

      do var = 0, nstatecomp - 1
         
      do i2 = iblo2,ibhi2
      do i1 = iblo1,ibhi1
      do i0 = iblo0,ibhi0

          slope (i0,i1,i2,var) = half * (
     &        state (i0+ii0,i1+ii1,i2+ii2,var) -
     &        state (i0-ii0,i1-ii1,i2-ii2,var) )
          
      enddo
      enddo
      enddo
       end do
      return
      end
      subroutine INTERPHISIDESLOPE(
     &           slope
     &           ,islopelo0,islopelo1,islopelo2
     &           ,islopehi0,islopehi1,islopehi2
     &           ,nslopecomp
     &           ,state
     &           ,istatelo0,istatelo1,istatelo2
     &           ,istatehi0,istatehi1,istatehi2
     &           ,nstatecomp
     &           ,iblo0,iblo1,iblo2
     &           ,ibhi0,ibhi1,ibhi2
     &           ,dir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nslopecomp
      integer islopelo0,islopelo1,islopelo2
      integer islopehi0,islopehi1,islopehi2
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1,
     &           islopelo2:islopehi2,
     &           0:nslopecomp-1)
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
      integer dir
      integer i0,i1,i2
      integer ii0,ii1,ii2
      integer var
      
      ii0=CHF_ID(0, dir)

      ii1=CHF_ID(1, dir)

      ii2=CHF_ID(2, dir)

      do var = 0, nstatecomp - 1
         
      do i2 = iblo2,ibhi2
      do i1 = iblo1,ibhi1
      do i0 = iblo0,ibhi0

          slope (i0,i1,i2,var) =
     &          state ( i0+ii0,i1+ii1,i2+ii2, var)
     &        - state ( i0,i1,i2, var)
          
      enddo
      enddo
      enddo
       enddo
      return
      end
      subroutine INTERPLOSIDESLOPE(
     &           slope
     &           ,islopelo0,islopelo1,islopelo2
     &           ,islopehi0,islopehi1,islopehi2
     &           ,nslopecomp
     &           ,state
     &           ,istatelo0,istatelo1,istatelo2
     &           ,istatehi0,istatehi1,istatehi2
     &           ,nstatecomp
     &           ,iblo0,iblo1,iblo2
     &           ,ibhi0,ibhi1,ibhi2
     &           ,dir
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nslopecomp
      integer islopelo0,islopelo1,islopelo2
      integer islopehi0,islopehi1,islopehi2
      REAL_T slope(
     &           islopelo0:islopehi0,
     &           islopelo1:islopehi1,
     &           islopelo2:islopehi2,
     &           0:nslopecomp-1)
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
      integer dir
      integer i0,i1,i2, var
      integer ii0,ii1,ii2
      
      ii0=CHF_ID(0, dir)

      ii1=CHF_ID(1, dir)

      ii2=CHF_ID(2, dir)

      do var = 0, nstatecomp - 1
         
      do i2 = iblo2,ibhi2
      do i1 = iblo1,ibhi1
      do i0 = iblo0,ibhi0

         slope (i0,i1,i2,var) =
     &        state (  i 0, i 1, i 2, var) -
     &        state (  i0-ii0, i1-ii1, i2-ii2, var)
          
      enddo
      enddo
      enddo
       end do
      return
      end
      subroutine INTERPLIMIT(
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
     &           ,state
     &           ,istatelo0,istatelo1,istatelo2
     &           ,istatehi0,istatehi1,istatehi2
     &           ,nstatecomp
     &           ,ibcoarselo0,ibcoarselo1,ibcoarselo2
     &           ,ibcoarsehi0,ibcoarsehi1,ibcoarsehi2
     &           ,ibnlo0,ibnlo1,ibnlo2
     &           ,ibnhi0,ibnhi1,ibnhi2
     &           ,iphysdomainlo0,iphysdomainlo1,iphysdomainlo2
     &           ,iphysdomainhi0,iphysdomainhi1,iphysdomainhi2
     &           )

      implicit none
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
      integer nstatecomp
      integer istatelo0,istatelo1,istatelo2
      integer istatehi0,istatehi1,istatehi2
      REAL_T state(
     &           istatelo0:istatehi0,
     &           istatelo1:istatehi1,
     &           istatelo2:istatehi2,
     &           0:nstatecomp-1)
      integer ibcoarselo0,ibcoarselo1,ibcoarselo2
      integer ibcoarsehi0,ibcoarsehi1,ibcoarsehi2
      integer ibnlo0,ibnlo1,ibnlo2
      integer ibnhi0,ibnhi1,ibnhi2
      integer iphysdomainlo0,iphysdomainlo1,iphysdomainlo2
      integer iphysdomainhi0,iphysdomainhi1,iphysdomainhi2
      integer  i0, i1, i2, var
      integer  ii0, ii1, ii2
      integer  in0, in1, in2
      REAL_T statemax, statemin, deltasum,  eta
      do var = 0, nislopecomp - 1
         
      do i2 = ibcoarselo2,ibcoarsehi2
      do i1 = ibcoarselo1,ibcoarsehi1
      do i0 = ibcoarselo0,ibcoarsehi0

             statemax = state ( i0,i1,i2, var )
             statemin = state ( i0,i1,i2, var )
             
      do ii2 = ibnlo2,ibnhi2
      do ii1 = ibnlo1,ibnhi1
      do ii0 = ibnlo0,ibnhi0

             
                 in0 = i0 + ii0
                 in1 = i1 + ii1
                 in2 = i2 + ii2
                 if (
                 
     &                in0 .ge. iphysdomainlo0 .and.
     &                in0 .le. iphysdomainhi0 
     &                .and.
     &                in1 .ge. iphysdomainlo1 .and.
     &                in1 .le. iphysdomainhi1 
     &                .and.
     &                in2 .ge. iphysdomainlo2 .and.
     &                in2 .le. iphysdomainhi2 
     &                ) then
                    statemax = max ( statemax, state(in0,in1,in2,var))
                    statemin = min ( statemin, state(in0,in1,in2,var))
                 endif
             
      enddo
      enddo
      enddo
             deltasum = half * (
                
     &            abs ( islope ( i0,i1,i2, var ) )
     &            +
     &            abs ( jslope ( i0,i1,i2, var ) )
     &            +
     &            abs ( kslope ( i0,i1,i2, var ) )
     &            )
#if CH_SPACEDIM > 3
                call MAYDAY_ERROR()
#endif
              eta = min(statemax - state(i0,i1,i2,var),
     &                  state(i0,i1,i2,var) - statemin)
              if( eta .le. 1.e-9*abs(statemax) ) then
                 eta = zero
              else
              if (deltasum .gt. eta) then
                eta = eta/deltasum
              else
                eta = one
              endif
              endif
              
              islope ( i0,i1,i2, var ) =
     &             eta * islope ( i0,i1,i2, var ) 
              jslope ( i0,i1,i2, var ) =
     &             eta * jslope ( i0,i1,i2, var ) 
              kslope ( i0,i1,i2, var ) =
     &             eta * kslope ( i0,i1,i2, var ) 
#if CH_SPACEDIM > 3
              call MAYDAY_ERROR()
#endif
         
      enddo
      enddo
      enddo
      end do
      return
      end
      subroutine INTERPLINEAR(
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
     &           ,dx
     &           ,stretch
     &           ,xbeg
     &           ,ixbeghi0
     &           ,ibreflo0,ibreflo1,ibreflo2
     &           ,ibrefhi0,ibrefhi1,ibrefhi2
     &           ,geometry
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
      REAL_T dx
      REAL_T stretch
      integer ixbeghi0
      REAL_T xbeg(
     &           0:ixbeghi0)
      integer ibreflo0,ibreflo1,ibreflo2
      integer ibrefhi0,ibrefhi1,ibrefhi2
      integer geometry
      integer ic0,ic1,ic2
      integer if0,if1,if2
      integer ii0,ii1,ii2
      integer var, id
      REAL_T  xlf0,xlf1,xlf2
      REAL_T  xrf0,xrf1,xrf2
      REAL_T  xlc0,xlc1,xlc2
      REAL_T  xrc0,xrc1,xrc2
      REAL_T dxf, volume
      if (geometry .eq. 1) then 
      do var = 0, nfinecomp - 1
          
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

              
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

              
                  if0 = ic0*ref_ratio + ii0
                  if1 = ic1*ref_ratio + ii1
                  if2 = ic2*ref_ratio + ii2
              
                  if (dir .eq. 0) then
                      id = ii0
                  else if (dir .eq. 1) then
                      id = ii1
                  else if (dir .eq. 2) then
                      id = ii2
                  endif
                  dxf = -half + ( (id+half) / ref_ratio )
                  fine( if0,if1,if2,var) =
     &                 fine( if0,if1,if2,var) +
     &                 dxf * slope (   ic 0,  ic 1,  ic 2, var )
              
      enddo
      enddo
      enddo
          
      enddo
      enddo
      enddo
      end do
      end if
      if ((geometry .eq. 2) .or. (geometry .eq. 5)) then
      do var = 0, nfinecomp - 1
          
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

               xlc0 = xbeg(0)+ic0*ref_ratio*dx
               xrc0 = xlc0 + ref_ratio*dx;
              
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

              
                  if0 = ic0*ref_ratio + ii0
                  volume = abs(xbeg(0)/dx+if0+half)
                  if1 = ic1*ref_ratio + ii1
                  if2 = ic2*ref_ratio + ii2
              
                  if (dir .eq. 0) then
                      xlf0 = xbeg(0)+if0*dx
                      xrf0 = xlf0 + dx
                      dxf = half*(xlf0*xlf0+xrf0*xrf0-xlc0*xlc0-xrc0*xrc0)/
     &                           (xrc0*xrc0-xlc0*xlc0)
                  else if (dir .eq. 1) then
                      dxf = -half + ( (ii1+half) / ref_ratio )
                  else if (dir .eq. 2) then
                      dxf = -half + ( (ii2+half) / ref_ratio )
                  endif
                  fine( if0,if1,if2,var) =
     &                 fine( if0,if1,if2,var) +
     &                 dxf * slope (   ic 0,  ic 1,  ic 2, var )*volume
              
      enddo
      enddo
      enddo
          
      enddo
      enddo
      enddo
      end do
      endif
      if (geometry .eq. 3) then
      do var = 0, nfinecomp - 1
          
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

              
               xlc0 = xbeg(0)+ic0*ref_ratio*dx
               xrc0 = xlc0 + ref_ratio*dx
               xlc1 = xbeg(1)+ic1*ref_ratio*dx*stretch
               xrc1 = xlc1 + ref_ratio*dx*stretch
                                 
              
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

              
                  if0 = ic0*ref_ratio + ii0
                  xlf0 = xbeg(0)+if0*dx
                  xrf0 = xlf0 + dx
                  volume = (xrf0*xrf0*xrf0-xlf0*xlf0*xlf0)*third
                  if1 = ic1*ref_ratio + ii1
                  xlf1 = xbeg(1)+if1*dx*stretch
                  xrf1 = xlf1 + dx*stretch
                  volume = volume*(cos(xlf1)-cos(xrf1))
                  if2 = ic2*ref_ratio + ii2
              
                  if (dir .eq. 0) then
                      dxf = half*(xlf0*xlf0*xlf0+xrf0*xrf0*xrf0-xlc0*xlc0*xlc0-xrc0*xrc0*xrc0)/
     &                           (xrc0*xrc0*xrc0-xlc0*xlc0*xlc0)
                  else if (dir .eq. 1) then
                      dxf = half*(cos(xlf1)+cos(xrf1)-cos(xlc1)-cos(xrc1))/
     &                           (cos(xrc1)-cos(xlc1))
                  else if (dir .eq. 2) then
                      dxf = -half + ( (ii2+half) / ref_ratio )
                  endif
                  fine( if0,if1,if2,var) =
     &                 fine( if0,if1,if2,var) +
     &                 dxf * slope (   ic 0,  ic 1,  ic 2, var )*volume
              
      enddo
      enddo
      enddo
          
      enddo
      enddo
      enddo
      end do
      endif
      if (geometry .eq. 4) then
      do var = 0, nfinecomp - 1
          
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

              
               xlc0 = ic0*ref_ratio*dx
               xrc0 = xlc0 + ref_ratio*dx
               xlc1 = xbeg(1)+ic1*ref_ratio*dx*stretch
               xrc1 = xlc1 + ref_ratio*dx*stretch
                                 
              
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

              
                  if0 = ic0*ref_ratio + ii0
                  xlf0 = if0*dx
                  xrf0 = xlf0 + dx
                  volume = (exp(three*xrf0)-exp(three*xlf0))*third
                  if1 = ic1*ref_ratio + ii1
                  xlf1 = xbeg(1)+if1*dx*stretch
                  xrf1 = xlf1 + dx*stretch
                  volume = volume*(cos(xlf1)-cos(xrf1))
                  if2 = ic2*ref_ratio + ii2
              
                  if (dir .eq. 0) then
                      dxf = half*(exp(three*xlf0)+exp(three*xrf0)-exp(three*xlc0)-exp(three*xrc0))/
     &                           (exp(three*xrc0)-exp(three*xlc0))
                  else if (dir .eq. 1) then
                      dxf = half*(cos(xlf1)+cos(xrf1)-cos(xlc1)-cos(xrc1))/
     &                           (cos(xrc1)-cos(xlc1))
                  else if (dir .eq. 2) then
                      dxf = -half + ( (ii2+half) / ref_ratio )
                  endif
                  fine( if0,if1,if2,var) =
     &                 fine( if0,if1,if2,var) +
     &                 dxf * slope (   ic 0,  ic 1,  ic 2, var )*volume
              
      enddo
      enddo
      enddo
          
      enddo
      enddo
      enddo
      end do
      endif
      if (geometry .eq. 6) then
      do var = 0, nfinecomp - 1
          
      do ic2 = iblo2,ibhi2
      do ic1 = iblo1,ibhi1
      do ic0 = iblo0,ibhi0

               xlc0 = ic0*ref_ratio*dx
               xrc0 = xlc0 + ref_ratio*dx;
              
      do ii2 = ibreflo2,ibrefhi2
      do ii1 = ibreflo1,ibrefhi1
      do ii0 = ibreflo0,ibrefhi0

              
                  if0 = ic0*ref_ratio + ii0
                  xlf0 = if0*dx
                  xrf0 = xlf0 + dx
                  volume = (exp(two*xrf0)-exp(two*xlf0))*half
                  if1 = ic1*ref_ratio + ii1
                  if2 = ic2*ref_ratio + ii2
              
                  if (dir .eq. 0) then
                      dxf = half*(exp(two*xlf0)+exp(two*xrf0)-exp(two*xlc0)-exp(two*xrc0))/
     &                           (exp(two*xrc0)-exp(two*xlc0))
                  else if (dir .eq. 1) then
                      dxf = -half + ( (ii1+half) / ref_ratio )
                  else if (dir .eq. 2) then
                      dxf = -half + ( (ii2+half) / ref_ratio )
                  endif
                  fine( if0,if1,if2,var) =
     &                 fine( if0,if1,if2,var) +
     &                 dxf * slope (   ic 0,  ic 1,  ic 2, var )*volume
              
      enddo
      enddo
      enddo
          
      enddo
      enddo
      enddo
      end do
      endif
      return
      end
      subroutine INTERPHOMO_OLD(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,x1
     &           ,dxCrse
     &           ,idir
     &           ,ihilo
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T x1
      REAL_T dxCrse
      integer idir
      integer ihilo
      REAL_T x2, denom, idenom, x, xsquared, m1, m2
      REAL_T q1, q2
      REAL_T pa, pb, a, b
      INTEGER ncomp,  n
      INTEGER ii0,ii1,ii2
      INTEGER i0,i1,i2
      x2 = half*(three*x1+dxCrse)
      denom = one-((x1+x2)/x1)
      idenom = one/(denom)
      x = two*x1
      xsquared = x*x
      m1 = one/(x1*x1)
      m2 = one/(x1*(x1-x2))
      q1 = one/(x1-x2)
      q2 = x1+x2
      ihilo = ihilo*(-1)
      ncomp = nphicomp
      
      ii0= ihilo*CHF_ID(0, idir)

      ii1= ihilo*CHF_ID(1, idir)

      ii2= ihilo*CHF_ID(2, idir)

      do n = 0, ncomp-1
          
      do i2 = iregionlo2,iregionhi2
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

          pa=phi(i0+2*ii0,i1+2*ii1,i2+2*ii2,n)
          pb=phi(i0+ii0,i1+ii1,i2+ii2,n)
          a=((pb-pa)*m1 - (pb)*m2)*idenom
          b=(pb)*q1 - a*q2
          phi(i0,i1,i2,n) = a*xsquared + b*x + pa
          
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine INTERPHOMO(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,x1
     &           ,dxCrse
     &           ,idir
     &           ,ihilo
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T x1
      REAL_T dxCrse
      integer idir
      integer ihilo
      REAL_T c1, c2
      REAL_T pa, pb
      INTEGER ncomp,  n
      INTEGER ii0,ii1,ii2
      INTEGER i0,i1,i2
      c1=two*(dxCrse-x1)/(dxCrse+x1)
      c2=   -(dxCrse-x1)/(dxCrse+three*x1)
      ihilo = ihilo*(-1)
      ncomp = nphicomp
      
      ii0= ihilo*CHF_ID(0, idir)

      ii1= ihilo*CHF_ID(1, idir)

      ii2= ihilo*CHF_ID(2, idir)

      do n = 0, ncomp-1
          
      do i2 = iregionlo2,iregionhi2
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

          pa=phi(i0+ii0,i1+ii1,i2+ii2,n)
          pb=phi(i0+2*ii0,i1+2*ii1,i2+2*ii2,n)
          phi(i0,i1,i2,n) = c1*pa + c2*pb
          
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine INTERPHOMOLINEAR(
     &           phi
     &           ,iphilo0,iphilo1,iphilo2
     &           ,iphihi0,iphihi1,iphihi2
     &           ,nphicomp
     &           ,iregionlo0,iregionlo1,iregionlo2
     &           ,iregionhi0,iregionhi1,iregionhi2
     &           ,x1
     &           ,dxCrse
     &           ,idir
     &           ,ihilo
     &           )

      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0,0,0,0,1,0 ,0,0,0,0,0,1 /


      integer nphicomp
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL_T phi(
     &           iphilo0:iphihi0,
     &           iphilo1:iphihi1,
     &           iphilo2:iphihi2,
     &           0:nphicomp-1)
      integer iregionlo0,iregionlo1,iregionlo2
      integer iregionhi0,iregionhi1,iregionhi2
      REAL_T x1
      REAL_T dxCrse
      integer idir
      integer ihilo
      INTEGER ncomp,  n
      INTEGER ii0,ii1,ii2
      INTEGER i0,i1,i2
      REAL_T pa, factor
      ihilo = ihilo*(-1)
      ncomp = nphicomp
      
      ii0= ihilo*CHF_ID(0, idir)

      ii1= ihilo*CHF_ID(1, idir)

      ii2= ihilo*CHF_ID(2, idir)

      factor = one - two*x1/(x1+dxCrse)
      do n = 0, ncomp-1
          
      do i2 = iregionlo2,iregionhi2
      do i1 = iregionlo1,iregionhi1
      do i0 = iregionlo0,iregionhi0

          pa=phi(i0+ii0,i1+ii1,i2+ii2,n)
          phi(i0,i1,i2,n) = factor*pa
          
      enddo
      enddo
      enddo
      enddo
      return
      end
