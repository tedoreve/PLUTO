MODULE Mod_Sync
!############################################################################
!#  This section contains the initialization subroutines for the synchrotron object
!#  as well as the subroutines related to the object itself.
!############################################################################

!### Modules to be Included###
USE Mod_Info
	
IMPLICIT NONE

PRIVATE :: F

!### Defines Synchrotron Type ###
TYPE :: Sync

	!### Array to store the emissivity ###
	REAL*8, ALLOCATABLE, DIMENSION(:) :: j

	!### Pointer to the emission frequencies ###
	REAL*8, POINTER :: nu(:)

	!### Useful Integers ###
	INTEGER :: nphbins

	!### Flag to see if synchrotron is on ###
	LOGICAL :: synchrotron=.FALSE.

END TYPE

CONTAINS

!### Initializes the Synchrotron Object ###
SUBROUTINE Sync_Init(synco,nu,nphbins)
	TYPE(Sync), INTENT(INOUT) :: synco

	REAL*8, TARGET :: nu(nphbins)
	INTEGER, INTENT(IN) :: nphbins

	!### Allocate space for the emissivity ###
	ALLOCATE(synco%j(nphbins))

	!### Store number of emission bins ###
	synco%nphbins=nphbins

	!### Set up pointer to frequency ###
	synco%nu => nu

END SUBROUTINE Sync_Init

!### Clean up Synchrotron Object ###
SUBROUTINE Sync_Destroy(synco)
	TYPE (Sync), INTENT(INOUT) :: synco

	!### Clean-up pointers ###
	NULLIFY(synco%nu)

	!### De-Allocate Arrays ###
	DEALLOCATE(synco%j)

END SUBROUTINE Sync_Destroy

!### Does the synchrotron calculation ###
!### The formalism in "Cosmic Ray Astrophysics" by R. Schlickeiser (2002) is used for this calculation (S02 from here on out). ###
!### To remove the /str of the synchrotron emissivity one needs to multiply by 8*\pi/3 for theta=pi/2 case, as the emission is non-isotropic. ###
!### Psuedo units means that those are the units if the constants were included.  The constants to be included are found at the end of the problem. ###
SUBROUTINE Sync_Calc(synco,cro,physo,ne,B,theta)
	TYPE(Sync), INTENT(INOUT) :: synco
	TYPE(CR), INTENT(IN) :: cro
	TYPE(Phys), INTENT(IN) :: physo

	REAL*8, INTENT(IN) :: ne, B, theta
	REAL*8 :: nup, nu0, istheta, stheta, xit, inu, A, D, jl, ju, slope
	INTEGER :: i, j

	!### Going to integrate it anyways so might as well zero it now. ###
	synco%j=0d0

	!### We don't need to run this routine if B equals 0 or sin(theta) .EQ. 0d0 ###
	if (B .EQ. 0d0 .OR. dsin(theta) .EQ. 0d0) then
		RETURN
	endif

	!### Plasma Frequency ###
	!### This allows for the calculation to damp the synchrotron spectrum below the Plasma Frequency ###
	nup=0.5d0/physo%pi*dsqrt(4d0*physo%pi*ne*physo%q**2/physo%me)

	!### Non-relativistic Gyrofrequency for the CR species ###
	nu0=0.5d0*dabs(cro%q)*B/(physo%pi*cro%m*physo%c)

	!### Inverse of sin(theta) ###
	stheta=dsin(theta)
	istheta=1d0/stheta

	!### Marches through emission frequencies ###
	do j=1,synco%nphbins
		!### Inverse of emission frequency ###
		inu=1d0/synco%nu(j)

		!### Work variables that are used frequently ###
		A=2d0*synco%nu(j)/(3d0*nu0)
		D=(nup*inu)**2
		
		!### Calculate synchrotron power for first momentum for use later ###
		!### x*istheta, where x is defined in S02 (4.1.4) ###
		xit=A*cro%ig2(1)*(1d0+cro%g2(1)*D)**(1.5d0)*istheta

		!### Eq. 4.1.7 of S02 excluding constants multiplied by dN/dE for the first bin, erg/s/Hz/cm^3/str/erg ###
		jl=F(xit)*cro%n(1)

		!### Integrates CR spectrum ###
		do i=1,cro%nbins-1
			!### x*istheta, where x is defined in S02 (4.1.4) ###
			xit=A*cro%ig2(i+1)*(1d0+cro%g2(i+1)*D)**(1.5d0)*istheta
			
			!### Eq. 4.1.7 of S02 excluding constants multiplied by dNdE for the next bin up, psuedo-units erg/s/Hz/cm^3/str/erg ###
			ju=F(xit)*cro%n(i+1)

			!### Check to save time ###
			if (jl .EQ. 0d0 .AND. ju .EQ. 0d0) then
				jl=ju
				CYCLE
			endif

			!### We will assume that the function produced by the product of the synchrotron power ###
			!### with the CR spectrum is a peice-wise powerlaw for the sake of integration. ###

			!### Thus first we will calculate the powerlaw slope ###
			slope=(dlog(MAX(ju,1d-100))-dlog(MAX(jl,1d-100)))/(dlog(cro%e(i+1))-dlog(cro%e(i)))

			!### Now we integrate over cosmic ray energy using the assumed powerlaw ###
			!### We need to test the slope to see if it is close to -1. ###
			!### If it is close to -1 then an edge case occurs and we need to change the integration ###
			!### synco%j should be in psuedo-units of ergs/cm^3/s/str/Hz ###
			!### This is excluding the normalization constants at this point which will put it into proper units. ###
			if (slope .LT. -.999999d0 .AND. slope .GT. -1.000001d0) then
				!### We will now calculate the normalization of the powerlaw ###
				if (jl .NE. 0d0) then
					jl=jl/cro%e(i)**slope
				else
					jl=ju/cro%e(i+1)**slope
				endif

				synco%j(j)=jl*dlog(cro%e(i+1)/cro%e(i))+synco%j(j)
			else
				synco%j(j)=1d0/(slope+1d0)*(ju*cro%e(i+1)-jl*cro%e(i))+synco%j(j)
			endif

			!### Store the upper value for the next update. ###
			jl=ju
		enddo
	enddo

	!### Now include the constants from S02 4.1.7, now we are actually in ergs/cm^3/s/str/Hz ###
	synco%j=2d0*physo%pi*dsqrt(3d0)*physo%r0*cro%m*physo%c*nu0*stheta*synco%j
			

END SUBROUTINE Sync_Calc

!### S02 (4.1.6c) ###
REAL*8 FUNCTION F(t)
	REAL*8 t

	F=1.78d0*t**(0.3d0)*dexp(-t)

END FUNCTION F

END MODULE Mod_Sync
