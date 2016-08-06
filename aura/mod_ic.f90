MODULE Mod_IC
!############################################################################
!#  This section contains the initialization subroutines for the 
!#  inverse-Compton (IC) object as well as the subroutines related to the 
!#  object itself.
!############################################################################

!### Modules to be Included###
USE Mod_Info
	
IMPLICIT NONE

PRIVATE :: G

!### Defines IC Type ###
TYPE :: IC

	!### Array to store the emissivity and work array ###
	REAL*8, ALLOCATABLE, DIMENSION(:) :: j, ea, iea

	!### Pointer to the emission frequencies ###
	REAL*8, POINTER :: nu(:), nupha(:), npha(:)

	!### Useful Integers ###
	INTEGER :: nphbins, nphabins

	!### Flag to see if IC is on ###
	LOGICAL :: ic=.FALSE.

END TYPE

CONTAINS

!### Initializes the IC Object ###
SUBROUTINE IC_Init(ico,nu,nphbins,nupha,npha,nphabins)
	TYPE(IC), INTENT(INOUT) :: ico

	REAL*8, TARGET :: nu(nphbins), nupha(nphabins), npha(nphabins)
	INTEGER, INTENT(IN) :: nphbins, nphabins

	!### Sanity Check ###
	if (nphabins .LT. 2) then
		write(*,*) 'IC ERROR: Number bins for the ambient photon field must be at least 2'
		STOP
	endif

	!### Allocate space for the emissivity and work array ###
	ALLOCATE(ico%j(nphbins),ico%ea(nphabins), ico%iea(nphabins))

	!### Store number of emission bins and ambient photon bins ###
	ico%nphbins=nphbins
	ico%nphabins=nphabins

	!### Set up pointer to frequency, ambient photon field and its frequency ###
	ico%nu => nu
	ico%nupha => nupha
	ico%npha => npha

END SUBROUTINE IC_Init

!### Clean up IC Object ###
SUBROUTINE IC_Destroy(ico)
	TYPE (IC), INTENT(INOUT) :: ico

	!### Clean-up pointers ###
	NULLIFY(ico%nu,ico%nupha,ico%npha)

	!### De-Allocate Arrays ###
	DEALLOCATE(ico%j,ico%ea, ico%iea)

END SUBROUTINE IC_Destroy

!### Does the IC calculation ###
!### The formalism in "Cosmic Ray Astrophysics" by R. Schlickeiser (2002) is used for this calculation (S02 from here on out). ###
SUBROUTINE IC_Calc(ico,cro,physo)
	TYPE(IC), INTENT(INOUT) :: ico
	TYPE(CR), INTENT(IN) :: cro
	TYPE(Phys), INTENT(IN) :: physo

	REAL*8 :: jl, ju, pl, pu, q, ge, mc2, imc2, eph, slope, qmin
	INTEGER :: i, j, k

	!### Going to integrate it anyways so might as well zero it now ###
	ico%j=0d0

	!### If there are no ambient photons then nothing needs to be done.
	if (MAXVAL(ico%npha) .EQ. 0d0 .OR. cro%m .NE. physo%me) then
		RETURN
	endif

	!### Useful quantities ###
	mc2=cro%m*physo%c**2
	imc2=1d0/mc2

	!### Setup work array for inverse of ambient photon energy. ###
	ico%ea=physo%h*ico%nupha
	ico%iea=1d0/ico%ea

	!### Marches through emission frequencies ###
	do k=1,ico%nphbins
		!### Current photon energy ###
		eph=physo%h*ico%nu(k)

		!### We will calculate the power at the lowest CR momentum so that we can do the integral in one pass rather than two. ###
		!### First we need to integrate over the ambient photon field ###
		jl=0d0

		!### Lower lower limit for q in the IC process from Eq. (4.2.4) in S02. ###
		qmin=0.25d0*cro%ig2(1)

		!### For the sake of speed we will start with the first photon bin ###
		!### S02 Eq. 4.2.3 ###
		ge=4d0*ico%ea(1)*cro%g(1)*imc2
		q=eph/(ge*(cro%g(1)*mc2-eph))

		!### Eq. 4.2.1 in S02 multiplied by the ambient photon field ###
		pl=ico%npha(1)*0.75d0*physo%sigmat*ico%iea(1)*cro%ig2(1)*G(q,ge,qmin)

		!### Integrates Ambient Photon Spectrum ###
		do i=1,ico%nphabins-1
			!### S02 Eq. 4.2.3 ###
			ge=4d0*ico%ea(i+1)*cro%g(1)*imc2
			q=eph/(ge*(cro%g(1)*mc2-eph))

			!### Eq. 4.2.1 in S02 multiplied by the ambient photon field ###
			pu=ico%npha(i+1)*0.75d0*physo%sigmat*ico%iea(i+1)*cro%ig2(1)*G(q,ge,qmin)

			!### Check to save time ###
			if (pl .EQ. 0d0 .AND. pu .EQ. 0d0) then
				pl=pu
				CYCLE
			endif

			!### We will assume that the function produced by the product of the cross-section ###
			!### with the ambient photon field is a peice-wise powerlaw for the sake of integration. ###

			!### Thus first we will calculate the powerlaw slope ###
			slope=(dlog(MAX(pu,1d-100))-dlog(MAX(pl,1d-100)))/(dlog(ico%ea(i+1))-dlog(ico%ea(i)))

			!### Now we integrate over ambient photon energy using the assumed powerlaw ###
			!### We need to test the slope to see if it is close to -1. ###
			!### If it is close to -1 then an edge case occurs and we need to change the integration ###
			!### ico%j should be in units of ergs/cm^3/s/str/erg ###
			if (slope .LT. -.999999d0 .AND. slope .GT. -1.000001d0) then
				!### We will now calculate the normalization of the powerlaw ###
				if (pl .NE. 0d0) then
					pl=pl/ico%ea(i)**slope
				else
					pl=pu/ico%ea(i+1)**slope
				endif

				jl=pl*dlog(ico%ea(i+1)/ico%ea(i))+jl
			else
				jl=1d0/(slope+1d0)*(pu*ico%ea(i+1)-pl*ico%ea(i))+jl
			endif

			!### Store the upper value for the next update. ###
			pl=pu
		enddo

		!### Complete Eq. 4.2.5 in S02 times the CR density for lowest momentum bin ###
		!### This is in erg/s/cm^3/erg/erg ###
		jl=physo%c*eph*jl*cro%n(1)
		
		!### Marches through particle energies ###
		do j=1,cro%nbins-1
			!### First we need to integrate over the ambient photon field ###
			ju=0d0

			!### Lower lower limit for q in the IC process from Eq. (4.2.4) in S02. ###
			qmin=0.25d0*cro%ig2(j+1)

			!### For the sake of speed we will start with the first photon bin ###
			!### S02 Eq. 4.2.3 ###
			ge=4d0*ico%ea(1)*cro%g(j+1)*imc2
			q=eph/(ge*(cro%g(j+1)*mc2-eph))

			!### Eq. 4.2.1 in S02 multiplied by the ambient photon field ###
			pl=ico%npha(1)*0.75d0*physo%sigmat*ico%iea(1)*cro%ig2(j+1)*G(q,ge,qmin)

			!### Integrates the Ambient Photon Spectrum ###
			do i=1,ico%nphabins-1
				!### S02 Eq. 4.2.3 ###
				ge=4d0*ico%ea(i+1)*cro%g(j+1)*imc2
				q=eph/(ge*(cro%g(j+1)*mc2-eph))

				!### Eq. 4.2.1 in S02 multiplied by the ambient photon field ###
				pu=ico%npha(i+1)*0.75d0*physo%sigmat*ico%iea(i+1)*cro%ig2(j+1)*G(q,ge,qmin)

				!### Check to save time ###
				if (pl .EQ. 0d0 .AND. pu .EQ. 0d0) then
					pl=pu
					CYCLE
				endif

				!### We will assume that the function produced by the product of the cross-section ###
				!### with the ambient photon field is a peice-wise powerlaw for the sake of integration. ###

				!### Thus first we will calculate the powerlaw slope ###
				slope=(dlog(MAX(pu,1d-100))-dlog(MAX(pl,1d-100)))/(dlog(ico%ea(i+1))-dlog(ico%ea(i)))

				!### Now we integrate over ambient photon energy using the assumed powerlaw ###
				!### We need to test the slope to see if it is close to -1. ###
				!### If it is close to -1 then an edge case occurs and we need to change the integration ###
				!### ico%j should be in units of ergs/cm^3/s/Hz ###
				if (slope .LT. -.999999d0 .AND. slope .GT. -1.000001d0) then
					!### We will now calculate the normalization of the powerlaw ###
					if (pl .NE. 0d0) then
						pl=pl/ico%ea(i)**slope
					else
						pl=pu/ico%ea(i+1)**slope
					endif

					ju=pl*dlog(ico%ea(i+1)/ico%ea(i))+ju
				else
					ju=1d0/(slope+1d0)*(pu*ico%ea(i+1)-pl*ico%ea(i))+ju
				endif

				!### Store the upper value for the next update. ###
				pl=pu
			enddo

			!### Complete Eq. 4.2.5 in S02 times the CR density ###
			!### This is in erg/s/cm^3/erg/erg ###
			ju=physo%c*eph*ju*cro%n(j+1)

			!### Check to save time ###
			if (jl .EQ. 0d0 .AND. ju .EQ. 0d0) then
				jl=ju
				CYCLE
			endif

			!### We will assume that the function produced by the product of the synchrotron power ###
			!### with the CR spectrum is a peice-wise powerlaw for the sake of integration. ###

			!### Thus first we will calculate the powerlaw slope ###
			slope=(dlog(MAX(ju,1d-100))-dlog(MAX(jl,1d-100)))/(dlog(cro%e(j+1))-dlog(cro%e(j)))

			!### Now we integrate over cosmic ray energy using the assumed powerlaw ###
			!### We need to test the slope to see if it is close to -1. ###
			!### If it is close to -1 then an edge case occurs and we need to change the integration ###
			!### ico%j should be in units of ergs/cm^3/s/erg ###
			if (slope .LT. -.999999d0 .AND. slope .GT. -1.000001d0) then
				!### We will now calculate the normalization of the powerlaw ###
				if (jl .NE. 0d0) then
					jl=jl/cro%e(j)**slope
				else
					jl=ju/cro%e(j+1)**slope
				endif

				ico%j(k)=jl*dlog(cro%e(j+1)/cro%e(j))+ico%j(k)
			else
				ico%j(k)=1d0/(slope+1d0)*(ju*cro%e(j+1)-jl*cro%e(j))+ico%j(k)
			endif

			!### Store the upper value for the next update. ###
			jl=ju
		enddo
	enddo

	!### S02 4.2.5 assumes the emission is isotropic and thus already folds in the factor 4*pi to take care of the angular component ###
	!### However in this code we will be working in units where the /str remains, so we will be dividing by 4*pi ###
	!### Also we need to convert from /erg to /Hz to get the proper emissivity units ###
	ico%j=0.25d0/physo%pi*physo%h*ico%j
	

END SUBROUTINE IC_Calc

!### Eq. 4.2.2 in S02 ###
REAL*8 FUNCTION G(q,ge,qmin)
	REAL*8 q, ge, qmin

	!### Evaluates bounds from Eq. 4.2.4 in S02.  Outside of this range IC is kinematically impossible so the cross-section is 0 ###
	if (q .GE. qmin .AND. q .LE. 1d0) then
		G=2d0*q*dlog(q)+(1d0+2d0*q)*(1d0-q)+0.5d0*(ge*q)**2*(1d0-q)/(1d0+ge*q)
	else
		G=0d0
	endif

END FUNCTION

END MODULE Mod_IC
