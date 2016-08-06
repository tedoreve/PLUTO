MODULE Mod_Brem
!############################################################################
!#  This section contains the initialization subroutines for the bremsstrahlung object
!#  as well as the subroutines related to the object itself.
!############################################################################

!### Modules to be Included###
USE Mod_Info
	
IMPLICIT NONE

PRIVATE :: Brem_phiH, Brem_phiHe

!### Defines Bremsstrahlung Type ###
TYPE :: Brem

	!### Array to store the emissivity ###
	REAL*8, ALLOCATABLE, DIMENSION(:) :: j

	!### Pointer to the emission frequencies ###
	REAL*8, POINTER :: nu(:)

	!### Useful Integers ###
	INTEGER :: nphbins

	!### Flag to see if bremsstrahlung is on ###
	LOGICAL :: bremsstrahlung=.FALSE.

END TYPE

CONTAINS

!### Initializes the Bremsstrahlung Object ###
SUBROUTINE Brem_Init(bremo,nu,nphbins)
	TYPE(Brem), INTENT(INOUT) :: bremo

	REAL*8, TARGET :: nu(nphbins)
	INTEGER, INTENT(IN) :: nphbins

	!### Allocate space for the emissivity ###
	ALLOCATE(bremo%j(nphbins))

	!### Store number of emission bins ###
	bremo%nphbins=nphbins

	!### Set up pointer to frequency ###
	bremo%nu => nu

END SUBROUTINE Brem_Init

!### Clean up Bremsstrahlung Object ###
SUBROUTINE Brem_Destroy(bremo)
	TYPE (Brem), INTENT(INOUT) :: bremo

	!### Clean-up pointers ###
	NULLIFY(bremo%nu)

	!### De-Allocate Arrays ###
	DEALLOCATE(bremo%j)

END SUBROUTINE Brem_Destroy

!### Does the Bremsstrahlung calculation ###
!### The formalism in "Cosmic Ray Astrophysics" by R. Schlickeiser (2002) is used the relativistic portion of this calculation (S02 from here on out). ###
!### Additionally for the non-relativistic end (Lorentz factor less than 2) the formalism in ###
!### "Classical Electrodynamics 3rd Ed." by J. D. Jackson 1999 was used for the cross-section (J99 from here on out) ###
!### Currently setup to handle an ambient medium containing ionized and neutral Hydrogen and Helium ###
SUBROUTINE Brem_Calc(bremo,cro,physo,ne,np,nalpha,nH,nHe)
	TYPE(Brem), INTENT(INOUT) :: bremo
	TYPE(CR), INTENT(IN) :: cro
	TYPE(Phys), INTENT(IN) :: physo

	REAL*8, INTENT(IN) :: ne, np, nalpha, nH, nHe
	REAL*8 :: jl, ju, slope, nion, eph, ieph, sigma, phiu, phi1, phi2, delta
	INTEGER :: i, j

	!### Going to integrate it anyways so might as well zero it now. ###
	bremo%j=0d0

	!### We will deal with the ionized protons and electrons first as they contribute to the bremsstrahlung in the same way. ###
	!### Total number of electron and protons
	nion=ne+np

	if (nion .NE. 0d0) then
		!### Marches through emission frequencies ###
		do j=1,bremo%nphbins
			!### Current photon energy and it's inverse ###
			eph=physo%h*bremo%nu(j)
			ieph=1d0/eph

			!### Calculate bremsstrahlung power for the first momentum for use later ###
			!### Bremsstrahlung cross section in units of cm^3 ###
			!### These are not the same units as S02 for the cross section ###
			if (cro%g(1) .LT. 2d0) then
				!### If the particle Lorentz factor is smaller than 2 we will consider it to be non-relativistic ###
				!### Hence we will use the form from J99 (15.29) which requires being divided by (4*pi)**2 ###
				if (cro%e(1) .GE. eph) then
					sigma=16d0/3d0*(2d0*physo%pi)/(physo%h*(4d0*physo%pi)**2)*physo%q**2/physo%c*physo%r0**2 &
					*1d0/(1d0-cro%ig2(1))*dlog((dsqrt(cro%e(1))+dsqrt(cro%e(1)-eph))**2*ieph)
				else
					sigma=0d0
				endif

				if (1d0/sigma .EQ. 0d0 .OR. ISNAN(sigma) .OR. (cro%e(1)-eph) .LT. 0d0 .OR. sigma .LT. 0d0) then
					sigma=0d0
				endif
			else
				!### Eq. 4.4.2 in S02 ###
				if (cro%te(1) .GE. eph) then
					phiu=4d0*(dlog(2d0*cro%g(1)*ieph*(cro%te(1)-eph))-0.5d0)
				else
					phiu=0d0
				endif

				!### Relativistic Cross section from S02 (4.4.1) ###
				sigma=3d0/physo%pi*0.125d0*physo%alpha*physo%sigmat*((1d0+(1d0-cro%ite(1)*eph)**2)*phiu-2d0/3d0*(1d0-eph*cro%ite(1))*phiu)
			endif

			!### The cross section cannot be less than 0 ###
			sigma=MAX(sigma,0d0)


			!### Eq. 4.4.4 in S02 multiplied by dN/dE for this ambient species erg/s/Hz/cm^3/erg ###
			!### Note that there appears to be a mistake in the units for S02.  What he lists is in ergs/s not ergs/s/eV. ###
			!### This includes the correction for if we are calculating bremsstrahlung for a CR species ###
			!### That is not electrons.  For a discussion see Dogiel et. al. A&A 382 (2002) 730 ###
			jl=physo%c*physo%h*nion*sigma*cro%n(1)*(cro%q/physo%q)**4*physo%me/cro%m

			!### Integrates the CR Spectrum ###
			do i=1,cro%nbins-1
				!### Bremsstrahlung cross section in units of cm^3 ###
				!### These are not the same units as S02 for the cross section ###
				if (cro%g(i+1) .LT. 2d0) then
					!### If the particle Lorentz factor is smaller than 2 we will consider it to be non-relativistic ###
					!### Hence we will use the form from J99 (15.29) which requires being divided by (4*pi)**2 ###
					if (cro%e(i+1) .GE. eph) then
						sigma=16d0/3d0*(2d0*physo%pi)/(physo%h*(4d0*physo%pi)**2)*physo%q**2/physo%c*physo%r0**2 &
						*1d0/(1d0-cro%ig2(i+1))*dlog((dsqrt(cro%e(i+1))+dsqrt(cro%e(i+1)-eph))**2*ieph)
					else
						sigma=0d0
					endif

					if (1d0/sigma .EQ. 0d0 .OR. ISNAN(sigma) .OR. (cro%e(i+1)-eph) .LT. 0d0 .OR. sigma .LT. 0d0) then
						sigma=0d0
					endif
				else
					!### Eq. 4.4.2 in S02 ###
					if (cro%te(i+1) .GE. eph) then
						phiu=4d0*(dlog(2d0*cro%g(i+1)*ieph*(cro%te(i+1)-eph))-0.5d0)
					else
						phiu=0d0
					endif

					!### Relativistic Cross section from S02 (4.4.1) ###
					sigma=3d0/physo%pi*0.125d0*physo%alpha*physo%sigmat*((1d0+(1d0-cro%ite(i+1)*eph)**2)*phiu-2d0/3d0*(1d0-eph*cro%ite(i+1))*phiu)

				endif

        	                !### The cross section cannot be less than 0 ###
	                        sigma=MAX(sigma,0d0)


				!### Eq. 4.4.4 in S02 multiplied by dN/dE for this species erg/s/Hz/cm^3/erg ###
				!### Note that there appears to be a mistake in the units for S02.  What he lists is in ergs/s not ergs/s/eV. ###
				!### This includes the correction for if we are calculating bremsstrahlung for a CR species ###
				!### That is not electrons.  For a discussion see Dogiel et. al. A&A 382 (2002) 730 ###
				ju=physo%c*physo%h*nion*sigma*cro%n(i+1)*(cro%q/physo%q)**4*physo%me/cro%m

				!### Check to save time ###
				if (jl .EQ. 0d0 .AND. ju .EQ. 0d0) then
					jl=ju
					CYCLE
				endif

				!### For this case it appears to work better if we use trapezoidal rule rather than a powerlaw assumption ###
				bremo%j(j)=0.5d0*(cro%e(i+1)-cro%e(i))*(jl+ju)+bremo%j(j)

				!### Store the upper value for the next update. ###
				jl=ju
			enddo
		enddo
	endif

	!### Next we will deal with alpha particles ###
	if (nalpha .NE. 0d0) then
		!### Marches through emission frequencies ###
		do j=1,bremo%nphbins
			!### Current photon energy and it's inverse ###
			eph=physo%h*bremo%nu(j)
			ieph=1d0/eph

			!### Calculate bremsstrahlung power for the first momentum for use later ###
			!### Bremsstrahlung cross section in units of cm^3 ###
			!### These are not the same units as S02 for the cross section ###
			if (cro%g(1) .LT. 2d0) then
				!### If the particle Lorentz factor is smaller than 2 we will consider it to be non-relativistic ###
				!### Hence we will use the form from J99 (15.29) which requires being divided by (4*pi)**2 ###
				if (cro%e(1) .GE. eph) then
					sigma=16d0*16d0/3d0*(2d0*physo%pi)/(physo%h*(4d0*physo%pi)**2)*physo%q**2/physo%c*physo%r0**2 &
					*1d0/(1d0-cro%ig2(1))*dlog((dsqrt(cro%e(1))+dsqrt(cro%e(1)-eph))**2*ieph)
				else
					sigma=0d0
				endif

				if (1d0/sigma .EQ. 0d0 .OR. ISNAN(sigma) .OR. (cro%e(1)-eph) .LT. 0d0 .OR. sigma .LT. 0d0) then
					sigma=0d0
				endif
			else
				!### Eq. 4.4.2 in S02 ###
				if (cro%te(1) .GE. eph) then
					phiu=16d0*4d0*(dlog(2d0*cro%g(1)*ieph*(cro%te(1)-eph))-0.5d0)
				else
					phiu=0d0
				endif

				!### Relativistic Cross section from S02 (4.4.1) ###
				sigma=3d0/physo%pi*0.125d0*physo%alpha*physo%sigmat*((1d0+(1d0-cro%ite(1)*eph)**2)*phiu-2d0/3d0*(1d0-eph*cro%ite(1))*phiu)
			endif

                        !### The cross section cannot be less than 0 ###
                        sigma=MAX(sigma,0d0)

			!### Eq. 4.4.4 in S02 multiplied by dN/dE for this ambient species erg/s/Hz/cm^3/erg ###
			!### This includes the correction for if we are calculating bremsstrahlung for a CR species ###
			!### Note that there appears to be a mistake in the units for S02.  What he lists is in ergs/s not ergs/s/eV. ###
			!### That is not electrons.  For a discussion see Dogiel et. al. A&A 382 (2002) 730 ###
			jl=physo%c*physo%h*nalpha*sigma*cro%n(1)*(cro%q/physo%q)**4*physo%me/cro%m

			!### Integrates the CR Spectrum ###
			do i=1,cro%nbins-1
				!### Bremsstrahlung cross section in units of cm^3 ###
				!### These are not the same units as S02 for the cross section ###
				if (cro%g(i+1) .LT. 2d0) then
					!### If the particle Lorentz factor is smaller than 2 we will consider it to be non-relativistic ###
					!### Hence we will use the form from J99 (15.29) which requires being divided by (4*pi)**2 ###
					if (cro%e(i+1) .GE. eph) then
						sigma=16d0*16d0/3d0*(2d0*physo%pi)/(physo%h*(4d0*physo%pi)**2)*physo%q**2/physo%c*physo%r0**2 &
						*1d0/(1d0-cro%ig2(i+1))*dlog((dsqrt(cro%e(i+1))+dsqrt(cro%e(i+1)-eph))**2*ieph)
					else
						sigma=0d0
					endif

					if (1d0/sigma .EQ. 0d0 .OR. ISNAN(sigma) .OR. (cro%e(i+1)-eph) .LT. 0d0 .OR. sigma .LT. 0d0) then
						sigma=0d0
					endif
				else
					!### Eq. 4.4.2 in S02 ###
					if (cro%te(i+1) .GE. eph) then
						phiu=16d0*4d0*(dlog(2d0*cro%g(i+1)*ieph*(cro%te(i+1)-eph))-0.5d0)
					else
						phiu=0d0
					endif

					!### Relativistic Cross section from S02 (4.4.1) ###
					sigma=3d0/physo%pi*0.125d0*physo%alpha*physo%sigmat*((1d0+(1d0-cro%ite(i+1)*eph)**2)*phiu-2d0/3d0*(1d0-eph*cro%ite(i+1))*phiu)
				endif

        	                !### The cross section cannot be less than 0 ###
	                        sigma=MAX(sigma,0d0)


				!### Eq. 4.4.4 in S02 multiplied by dN/dE for this species erg/s/Hz/cm^3/erg ###
				!### This includes the correction for if we are calculating bremsstrahlung for a CR species ###
				!### That is not electrons.  For a discussion see Dogiel et. al. A&A 382 (2002) 730 ###
				!### Note that there appears to be a mistake in the units for S02.  What he lists is in ergs/s not ergs/s/eV. ###
				ju=physo%c*physo%h*nalpha*sigma*cro%n(i+1)*(cro%q/physo%q)**4*physo%me/cro%m

				!### Check to save time ###
				if (jl .EQ. 0d0 .AND. ju .EQ. 0d0) then
					jl=ju
					CYCLE
				endif

				!### For this case it appears to work better if we use trapezoidal rule rather than a powerlaw assumption ###
				bremo%j(j)=0.5d0*(cro%e(i+1)-cro%e(i))*(jl+ju)+bremo%j(j)

				!### Store the upper value for the next update. ###
				jl=ju
			enddo
		enddo
	endif

	!### Now we will handle the neutral Hydrogen ###
	if (nH .NE. 0d0) then
		!### Marches through emission frequencies ###
		do j=1,bremo%nphbins
			!### Current photon energy and it's inverse ###
			eph=physo%h*bremo%nu(j)
			ieph=1d0/eph

			!### Calculate bremsstrahlung power for the first momentum for use later ###
			!### Eq. 4.4.3 in S02 ###
			delta=0.25d0*eph*cro%ig(1)/(physo%alpha*(cro%te(1)-eph))

			!### Calculate phi_1 and phi_2 based on Delta ###
			if (delta .GT. 2d0) then
				!### If delta is greater than 2 then phi1 and phi1 reduce to phiu ###
				!### Eq. 4.4.2 in S02 ###
				phi1=4d0*(dlog(2d0*cro%g(1)*ieph*(cro%te(1)-eph))-0.5d0)
				phi2=phi1
			elseif (delta .GE. 0d0 .AND. delta .LE. 2d0) then

				!### Calculates phi1 and phi2 ###
				CALL Brem_phiH(delta,phi1,phi2)

			else
				phi1=0d0
				phi2=0d0
			endif

			!### Bremsstrahlung cross section in units of cm^3 ###
			!### These are not the same units as S02 for the cross section ###
			!### Relativistic Cross section from S02 (4.4.1) ###
			sigma=3d0/physo%pi*0.125d0*physo%alpha*physo%sigmat*((1d0+(1d0-cro%ite(1)*eph)**2)*phi1-2d0/3d0*(1d0-eph*cro%ite(1))*phi2)

                        !### The cross section cannot be less than 0 ###
                        sigma=MAX(sigma,0d0)

			!### Eq. 4.4.4 in S02 multiplied by dN/dE for this ambient species erg/s/Hz/cm^3/erg ###
			!### This includes the correction for if we are calculating bremsstrahlung for a CR species ###
			!### That is not electrons.  For a discussion see Dogiel et. al. A&A 382 (2002) 730 ###
			!### Note that there appears to be a mistake in the units for S02.  What he lists is in ergs/s not ergs/s/eV. ###
			jl=physo%c*physo%h*nH*sigma*cro%n(1)*(cro%q/physo%q)**4*physo%me/cro%m

			!### Integrates the CR Spectrum ###
			do i=1,cro%nbins-1
				!### Eq. 4.4.3 in S02 ###
				delta=0.25d0*eph*cro%ig(i+1)/(physo%alpha*(cro%te(i+1)-eph))

				!### Calculate phi_1 and phi_2 based on Delta ###
				if (delta .GT. 2d0) then
					!### If delta is greater than 2 then phi1 and phi1 reduce to phiu ###
					!### Eq. 4.4.2 in S02 ###
					phi1=4d0*(dlog(2d0*cro%g(i+1)*ieph*(cro%te(i+1)-eph))-0.5d0)
					phi2=phi1
				elseif (delta .GE. 0d0 .AND. delta .LE. 2d0) then

					!### Calculates phi1 and phi2 ###
					CALL Brem_phiH(delta,phi1,phi2)

				else
					phi1=0d0
					phi2=0d0
				endif

				!### Bremsstrahlung cross section in units of cm^3 ###
				!### These are not the same units as S02 for the cross section ###
				!### Relativistic Cross section from S02 (4.4.1) ###
				sigma=3d0/physo%pi*0.125d0*physo%alpha*physo%sigmat*((1d0+(1d0-cro%ite(i+1)*eph)**2)*phi1-2d0/3d0*(1d0-eph*cro%ite(i+1))*phi2)

        	                !### The cross section cannot be less than 0 ###
	                        sigma=MAX(sigma,0d0)


				!### Eq. 4.4.4 in S02 multiplied by dN/dE for this species erg/s/Hz/cm^3/erg ###
				!### This includes the correction for if we are calculating bremsstrahlung for a CR species ###
				!### That is not electrons.  For a discussion see Dogiel et. al. A&A 382 (2002) 730 ###
				!### Note that there appears to be a mistake in the units for S02.  What he lists is in ergs/s not ergs/s/eV. ###
				ju=physo%c*physo%h*nH*sigma*cro%n(i+1)*(cro%q/physo%q)**4*physo%me/cro%m

				!### Check to save time ###
				if (jl .EQ. 0d0 .AND. ju .EQ. 0d0) then
					jl=ju
					CYCLE
				endif

				!### For this case it appears to work better if we use trapezoidal rule rather than a powerlaw assumption ###
				bremo%j(j)=0.5d0*(cro%e(i+1)-cro%e(i))*(jl+ju)+bremo%j(j)

				!### Store the upper value for the next update. ###
				jl=ju
			enddo
		enddo
	endif

	!### Now we will handle the neutral Helium ###
	if (nHe .NE. 0d0) then
		!### Marches through emission frequencies ###
		do j=1,bremo%nphbins
			!### Current photon energy and it's inverse ###
			eph=physo%h*bremo%nu(j)
			ieph=1d0/eph

			!### Calculate bremsstrahlung power for the first momentum for use later ###
			!### Eq. 4.4.3 in S02 ###
			delta=0.25d0*eph*cro%ig(1)/(physo%alpha*(cro%te(1)-eph))

			!### Calculate phi_1 and phi_2 based on Delta ###
			if (delta .GT. 2d0) then
				!### If delta is greater than 2 then phi1 and phi1 reduce to phiu ###
				!### Eq. 4.4.2 in S02 ###
				phi1=16d0*4d0*(dlog(2d0*cro%g(1)*ieph*(cro%te(1)-eph))-0.5d0)
				phi2=phi1
			elseif (delta .GE. 0d0 .AND. delta .LE. 2d0) then

				!### Calculates phi1 and phi2 ###
				CALL Brem_phiHe(delta,phi1,phi2)

			else
				phi1=0d0
				phi2=0d0
			endif

			!### Bremsstrahlung cross section in units of cm^3 ###
			!### These are not the same units as S02 for the cross section ###
			!### Relativistic Cross section from S02 (4.4.1) ###
			sigma=3d0/physo%pi*0.125d0*physo%alpha*physo%sigmat*((1d0+(1d0-cro%ite(1)*eph)**2)*phi1-2d0/3d0*(1d0-eph*cro%ite(1))*phi2)

                        !### The cross section cannot be less than 0 ###
                        sigma=MAX(sigma,0d0)

			!### Eq. 4.4.4 in S02 multiplied by dN/dE for this ambient species erg/s/Hz/cm^3/erg ###
			!### This includes the correction for if we are calculating bremsstrahlung for a CR species ###
			!### That is not electrons.  For a discussion see Dogiel et. al. A&A 382 (2002) 730 ###
			!### Note that there appears to be a mistake in the units for S02.  What he lists is in ergs/s not ergs/s/eV. ###
			jl=physo%c*physo%h*nHe*sigma*cro%n(1)*(cro%q/physo%q)**4*physo%me/cro%m

			!### Integrates the CR Spectrum ###
			do i=1,cro%nbins-1
				!### Eq. 4.4.3 in S02 ###
				delta=0.25d0*eph*cro%ig(i+1)/(physo%alpha*(cro%te(i+1)-eph))

				!### Calculate phi_1 and phi_2 based on Delta ###
				if (delta .GT. 2d0) then
					!### If delta is greater than 2 then phi1 and phi1 reduce to phiu ###
					!### Eq. 4.4.2 in S02 ###
					phi1=16d0*4d0*(dlog(2d0*cro%g(i+1)*ieph*(cro%te(i+1)-eph))-0.5d0)
					phi2=phi1
				elseif (delta .GE. 0d0 .AND. delta .LE. 2d0) then

					!### Calculates phi1 and phi2 ###
					CALL Brem_phiHe(delta,phi1,phi2)

				else
					phi1=0d0
					phi2=0d0
				endif

				!### Bremsstrahlung cross section in units of cm^3 ###
				!### These are not the same units as S02 for the cross section ###
				!### Relativistic Cross section from S02 (4.4.1) ###
				sigma=3d0/physo%pi*0.125d0*physo%alpha*physo%sigmat*((1d0+(1d0-cro%ite(i+1)*eph)**2)*phi1-2d0/3d0*(1d0-eph*cro%ite(i+1))*phi2)

        	                !### The cross section cannot be less than 0 ###
	                        sigma=MAX(sigma,0d0)

				!### Eq. 4.4.4 in S02 multiplied by dN/dE for this species erg/s/Hz/cm^3/erg ###
				!### This includes the correction for if we are calculating bremsstrahlung for a CR species ###
				!### That is not electrons.  For a discussion see Dogiel et. al. A&A 382 (2002) 730 ###
				!### Note that there appears to be a mistake in the units for S02.  What he lists is in ergs/s not ergs/s/eV. ###
				ju=physo%c*physo%h*nHe*sigma*cro%n(i+1)*(cro%q/physo%q)**4*physo%me/cro%m

				!### Check to save time ###
				if (jl .EQ. 0d0 .AND. ju .EQ. 0d0) then
					jl=ju
					CYCLE
				endif

				!### For this case it appears to work better if we use trapezoidal rule rather than a powerlaw assumption ###
				bremo%j(j)=0.5d0*(cro%e(i+1)-cro%e(i))*(jl+ju)+bremo%j(j)

				!### Store the upper value for the next update. ###
				jl=ju
			enddo
		enddo
	endif

	!### S02 4.4.4 assumes the emission is isotropic and thus already folds in the factor 4*pi to take care of the angular component ###
	!### However in this code we will be working in units where the /str remains, so we will be dividing by 4*pi ###
	bremo%j=0.25d0/physo%pi*bremo%j
	
END SUBROUTINE Brem_Calc

!### Calculates phi1 and phi2 for Hydrogen based on Table 4.1 in S02 ###
SUBROUTINE Brem_phiH(delta,phi1,phi2)
	REAL*8, INTENT(IN) :: delta
	REAL*8, INTENT(OUT) :: phi1, phi2
	REAL*8 :: d(9), p1(9), p2(9), slope
	INTEGER :: i

	!### Column 1, 2, 3 of Table 4.1 S02 ###
	d = (/0d0, 0.01d0, 0.02d0, 0.05d0, 0.1d0, 0.2d0, 0.5d0, 1d0, 2d0/)
	p1= (/45.79d0, 45.43d0, 45.09d0, 44.11d0, 42.64d0, 40.16d0, 34.97d0, 29.97d0, 24.73d0/)
	p2= (/44.46d0, 44.38d0, 44.24d0, 43.65d0, 42.49d0, 40.19d0, 34.93d0, 29.78d0, 24.34d0/) 

	!### Now we will figure out in between which two bins delta lies. ###
	do i=1,8
		if (delta .GE. d(i) .AND. delta .LE. d(i+1)) EXIT
	enddo

	!### We will use simple linear interpolation for the inter bin values as the function is slowly varying. ###
	!### phi_1 ###
	slope=(p1(i+1)-p1(i))/(d(i+1)-d(i))
	phi1=slope*(delta-d(i))+p1(i)

	!### phi_2 ###
	slope=(p2(i+1)-p2(i))/(d(i+1)-d(i))
	phi2=slope*(delta-d(i))+p2(i)

END SUBROUTINE Brem_phiH

!### Calculates phi1 and phi2 for Helium based on Table 4.1 in S02 ###
SUBROUTINE Brem_phiHe(delta,phi1,phi2)
	REAL*8, INTENT(IN) :: delta
	REAL*8, INTENT(OUT) :: phi1, phi2
	REAL*8 :: d(9), p1(9), p2(9), slope
	INTEGER :: i

	!### Column 1, 4, 5 of Table 4.1 S02 ###
	d = (/0d0, 0.01d0, 0.02d0, 0.05d0, 0.1d0, 0.2d0, 0.5d0, 1d0, 2d0/)
	p1= (/134.6d0, 133.85d0, 133.11d0, 130.86d0, 127.17d0, 120.35d0, 104.60d0, 89.94d0, 74.19d0/)
	p2= (/131.4d0, 130.51d0, 130.33d0, 129.26d0, 126.76d0, 120.80d0, 105.21d0, 89.46d0, 73.03d0/)

	!### Now we will figure out in between which two bins delta lies. ###
	do i=1,8
		if (delta .GE. d(i) .AND. delta .LE. d(i+1)) EXIT
	enddo

	!### We will use simple linear interpolation for the inter bin values as the function is slowly varying. ###
	!### phi_1 ###
	slope=(p1(i+1)-p1(i))/(d(i+1)-d(i))
	phi1=slope*(delta-d(i))+p1(i)

	!### phi_2 ###
	slope=(p2(i+1)-p2(i))/(d(i+1)-d(i))
	phi2=slope*(delta-d(i))+p2(i)

END SUBROUTINE Brem_phiHe

END MODULE Mod_Brem
