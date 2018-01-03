      subroutine REGCORRECTTANVEL(
     & vel
     & ,ivello0,ivello1,ivello2
     & ,ivelhi0,ivelhi1,ivelhi2
     & ,grad
     & ,igradlo0,igradlo1,igradlo2
     & ,igradhi0,igradhi1,igradhi2
     & ,iinteriorboxlo0,iinteriorboxlo1,iinteriorboxlo2
     & ,iinteriorboxhi0,iinteriorboxhi1,iinteriorboxhi2
     & ,veldir
     & ,graddir
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer ivello0,ivello1,ivello2
      integer ivelhi0,ivelhi1,ivelhi2
      REAL*8 vel(
     & ivello0:ivelhi0,
     & ivello1:ivelhi1,
     & ivello2:ivelhi2)
      integer igradlo0,igradlo1,igradlo2
      integer igradhi0,igradhi1,igradhi2
      REAL*8 grad(
     & igradlo0:igradhi0,
     & igradlo1:igradhi1,
     & igradlo2:igradhi2)
      integer iinteriorboxlo0,iinteriorboxlo1,iinteriorboxlo2
      integer iinteriorboxhi0,iinteriorboxhi1,iinteriorboxhi2
      integer veldir
      integer graddir
      integer ioffvel , joffvel , koffvel
      integer ioffgrad, joffgrad, koffgrad
      integer i,j,k
      REAL*8 correction, factor
      ioffvel = chf_id(0,veldir)
      joffvel = chf_id(1,veldir)
      koffvel = chf_id(2,veldir)
      ioffgrad = chf_id(0,graddir)
      joffgrad = chf_id(1,graddir)
      koffgrad = chf_id(2,graddir)
      factor = (1.0d0)/(4.0d0)
      do k = iinteriorboxlo2,iinteriorboxhi2
      do j = iinteriorboxlo1,iinteriorboxhi1
      do i = iinteriorboxlo0,iinteriorboxhi0
      correction = factor*
     $ ( grad(i ,j ,k )
     $ + grad(i+ioffgrad-ioffvel,j+joffgrad-joffvel,k+koffgrad-koffvel)
     $ + grad(i -ioffvel,j -joffvel,k -koffvel)
     $ + grad(i+ioffgrad ,j+joffgrad ,k+koffgrad ))
      vel(i,j,k) = vel(i,j,k) - correction
      enddo
      enddo
      enddo
      return
      end
      subroutine MACDIVERGEF(
     & idcalclo0,idcalclo1,idcalclo2
     & ,idcalchi0,idcalchi1,idcalchi2
     & ,divf
     & ,idivflo0,idivflo1,idivflo2
     & ,idivfhi0,idivfhi1,idivfhi2
     & ,ndivfcomp
     & ,flux
     & ,ifluxlo0,ifluxlo1,ifluxlo2
     & ,ifluxhi0,ifluxhi1,ifluxhi2
     & ,nfluxcomp
     & ,facedir
     & ,nconserved
     & ,dx
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer idcalclo0,idcalclo1,idcalclo2
      integer idcalchi0,idcalchi1,idcalchi2
      integer ndivfcomp
      integer idivflo0,idivflo1,idivflo2
      integer idivfhi0,idivfhi1,idivfhi2
      REAL*8 divf(
     & idivflo0:idivfhi0,
     & idivflo1:idivfhi1,
     & idivflo2:idivfhi2,
     & 0:ndivfcomp-1)
      integer nfluxcomp
      integer ifluxlo0,ifluxlo1,ifluxlo2
      integer ifluxhi0,ifluxhi1,ifluxhi2
      REAL*8 flux(
     & ifluxlo0:ifluxhi0,
     & ifluxlo1:ifluxhi1,
     & ifluxlo2:ifluxhi2,
     & 0:nfluxcomp-1)
      integer facedir
      integer nconserved
      REAL*8 dx
      integer i, j, k
      integer ioff, joff, koff
      integer spacedim,iv
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      spacedim = 3
      do iv = 0,nconserved - 1
      do k = idcalclo2,idcalchi2
      do j = idcalclo1,idcalchi1
      do i = idcalclo0,idcalchi0
         divf(i,j,k,iv) = divf(i,j,k,iv) +
     & (flux(i+ioff,j+joff,k+koff,iv)
     & -flux(i ,j ,k ,iv))/dx
      enddo
      enddo
      enddo
      enddo
      return
      end
      subroutine MACGRADPHI(
     & gradphi
     & ,igradphilo0,igradphilo1,igradphilo2
     & ,igradphihi0,igradphihi1,igradphihi2
     & ,phi
     & ,iphilo0,iphilo1,iphilo2
     & ,iphihi0,iphihi1,iphihi2
     & ,facedir
     & ,dx
     & ,ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
     & ,ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
     & )
      implicit none
      integer CHF_ID(0:5,0:5)
      data CHF_ID/ 1,0,0,0,0,0 ,0,1,0,0,0,0 ,0,0,1,0,0,0 ,0,0,0,1,0,0 ,0
     &,0,0,0,1,0 ,0,0,0,0,0,1 /
      integer igradphilo0,igradphilo1,igradphilo2
      integer igradphihi0,igradphihi1,igradphihi2
      REAL*8 gradphi(
     & igradphilo0:igradphihi0,
     & igradphilo1:igradphihi1,
     & igradphilo2:igradphihi2)
      integer iphilo0,iphilo1,iphilo2
      integer iphihi0,iphihi1,iphihi2
      REAL*8 phi(
     & iphilo0:iphihi0,
     & iphilo1:iphihi1,
     & iphilo2:iphihi2)
      integer facedir
      REAL*8 dx
      integer ifaceboxlo0,ifaceboxlo1,ifaceboxlo2
      integer ifaceboxhi0,ifaceboxhi1,ifaceboxhi2
      integer i,j,k
      integer ioff,joff,koff
      ioff = chf_id(0,facedir)
      joff = chf_id(1,facedir)
      koff = chf_id(2,facedir)
      do k = ifaceboxlo2,ifaceboxhi2
      do j = ifaceboxlo1,ifaceboxhi1
      do i = ifaceboxlo0,ifaceboxhi0
      gradphi(i,j,k) =
     & ( phi(i ,j ,k )
     & - phi(i-ioff,j-joff,k-koff)
     & )/dx
      enddo
      enddo
      enddo
      return
      end
