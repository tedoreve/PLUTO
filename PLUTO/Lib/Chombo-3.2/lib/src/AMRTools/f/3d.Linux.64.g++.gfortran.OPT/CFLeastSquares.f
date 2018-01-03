      REAL*8 function xpower1d0avg(lower, upper, p)
      implicit none
      REAL*8 lower, upper
      integer p
      REAL*8 avg
      avg = (upper**(p+1) - lower**(p+1)) / (upper - lower)
      if ((p .gt. 0) .and. (mod(p, 2) .eq. 0)) then
         avg = avg - ((0.500d0)**p)
      endif
      avg = avg / ((p + 1) * (1.0d0))
      xpower1d0avg = avg
      return
      end
      REAL*8 function xpower1dcoarseind0avg(off, p)
      implicit none
      integer off, p
      REAL*8 lower, upper, avg
      REAL*8 xpower1d0avg
      lower = (off*(1.0d0)) - (0.500d0)
      upper = (off*(1.0d0)) + (0.500d0)
      avg = xpower1d0avg(lower, upper, p)
      xpower1dcoarseind0avg = avg
      return
      end
      subroutine xpowercoarseind0avg(
     &     nnbrs, npwrs, nbrs, pwrs, xpower)
      implicit none
      integer nnbrs, npwrs
      integer nbrs(nnbrs, 3), pwrs(npwrs, 3)
      REAL*8 xpower(nnbrs, npwrs)
      integer inbr, ipwr, nbr, pwr, idim
      REAL*8 xpower1dcoarseind0avg
      do inbr = 1, nnbrs
         do ipwr = 1, npwrs
            xpower(inbr, ipwr) = (1.0d0)
         end do
      end do
      do idim = 1, 3
         do inbr = 1, nnbrs
            nbr = nbrs(inbr, idim)
            do ipwr = 1, npwrs
               pwr = pwrs(ipwr, idim)
               xpower(inbr, ipwr) = xpower(inbr, ipwr) *
     &              xpower1dcoarseind0avg(nbr, pwr)
            end do
         end do
      end do
      return
      end
      REAL*8 function xpower1dfineind0avg(off, p, nref)
      implicit none
      integer off, p, nref
      REAL*8 lower, upper, avg
      REAL*8 xpower1d0avg
      lower = ( off     *(1.0d0)) / (nref*(1.0d0)) - (0.500d0)
      upper = ((off + 1)*(1.0d0)) / (nref*(1.0d0)) - (0.500d0)
      avg = xpower1d0avg(lower, upper, p)
      xpower1dfineind0avg = avg
      return
      end
      subroutine xpowerfineind0avg(
     &     nfine, npwrs, nref, pwrs, fine, xpower)
      implicit none
      integer nfine, npwrs, nref
      integer pwrs(npwrs, 3)
      integer fine(nfine, 3)
      REAL*8 xpower(nfine, npwrs)
      integer ifine, ipwr, pwr, idim, off
      REAL*8 xpower1dfineind0avg
      do ifine = 1, nfine
         do ipwr = 1, npwrs
            xpower(ifine, ipwr) = (1.0d0)
         end do
      end do
      do idim = 1, 3
         do ifine = 1, nfine
            off = fine(ifine, idim)
            do ipwr = 1, npwrs
               pwr = pwrs(ipwr, idim)
               xpower(ifine, ipwr) = xpower(ifine, ipwr) *
     &              xpower1dfineind0avg(off, pwr, nref)
            end do
         end do
      end do
      return
      end
      integer function findnbr0(nnbrs, nbrs)
      implicit none
      integer nnbrs
      integer nbrs(nnbrs, 3)
      integer inbr, idim
      logical found
      do inbr = 1, nnbrs
         found = .true.
         do idim = 1, 3
            found = found .and. (nbrs(inbr, idim) .eq. 0)
         end do
         if (found) then
            findnbr0 = inbr
         endif
      end do
      return
      end
      subroutine coarsefineleastsquares(
     &     nref, nnbrs, npwrs, nfine,
     &     nbrs, pwrs, fine, mat)
      implicit none
      integer nref, nnbrs, npwrs, nfine
      integer nbrs(nnbrs, 3)
      integer pwrs(npwrs, 3)
      integer fine(nfine, 3)
      REAL*8 mat(nfine, nnbrs)
      REAL*8 xpowc(nnbrs, npwrs)
      REAL*8 ata(npwrs, npwrs)
      REAL*8 work(npwrs)
      REAL*8 coeffs(npwrs, nnbrs)
      REAL*8 sumrow(npwrs)
      REAL*8 xpowf(nfine, npwrs)
      integer ipiv(npwrs), info, nbr0
      integer inbr, ipwr, ifine
      character*1 prec
      integer findnbr0
      prec = 'D'
      call xpowercoarseind0avg(nnbrs, npwrs, nbrs, pwrs, xpowc)
      call DGEMM('T', 'N', npwrs, npwrs, nnbrs, (1.0d0), xpowc,
     &     nnbrs, xpowc, nnbrs, (0.0d0), ata, npwrs)
      call DGETRF(npwrs, npwrs, ata, npwrs, ipiv, info)
      if (info .lt. 0) then
         print *, 'Error:  LAPACK ', prec, 'GETRF returned ', info
         print *, 'Illegal value in argument number ', -info
         call MAYDAY_ERROR()
      elseif (info .gt. 0) then
         print *, 'Error:  LAPACK ', prec, 'GETRF returned ', info
         print *, 'Singular matrix, 0 at diagonal element ', info
         call MAYDAY_ERROR()
      endif
      call DGETRI(npwrs, ata, npwrs, ipiv, work, npwrs, info)
      if (info .lt. 0) then
         print *, 'Error:  LAPACK ', prec, 'GETRI returned ', info
         print *, 'Illegal value in argument number ', -info
         call MAYDAY_ERROR()
      elseif (info .gt. 0) then
         print *, 'Error:  LAPACK ', prec, 'GETRI returned ', info
         print *, 'Singular matrix, 0 at diagonal element ', info
         call MAYDAY_ERROR()
      endif
      call DGEMM('N', 'T', npwrs, nnbrs, npwrs, (1.0d0), ata,
     &     npwrs, xpowc, nnbrs,0, coeffs, npwrs)
      do ipwr = 1, npwrs
         sumrow(ipwr) = 0
         do inbr = 1, nnbrs
            sumrow(ipwr) = sumrow(ipwr) + coeffs(ipwr, inbr)
         end do
      end do
      nbr0 = findnbr0(nnbrs, nbrs)
      do ipwr = 1, npwrs
         coeffs(ipwr, nbr0) = coeffs(ipwr, nbr0) - sumrow(ipwr)
      end do
      call xpowerfineind0avg(nfine, npwrs, nref, pwrs, fine, xpowf)
      call DGEMM('N', 'N', nfine, nnbrs, npwrs, (1.0d0), xpowf,
     &     nfine, coeffs, npwrs, 0, mat, nfine)
      do ifine = 1, nfine
         mat(ifine, nbr0) = mat(ifine, nbr0) + (1.0d0)
      end do
      return
      end
