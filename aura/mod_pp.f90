MODULE Mod_PP
!############################################################################
!#  This section contains the initialization subroutines for the proton-proton interaction object
!#  as well as the subroutines related to the object itself.
!############################################################################

!### Modules to be Included ###
USE Mod_Info
	
IMPLICIT NONE

PRIVATE :: sigmai, Fg, Fe, Fnum

!### Defines proton-proton interaction Type ###
TYPE :: PP

	TYPE(CR) :: eo !Secondary Electrons
	TYPE(CR) :: po !Secondary Positrons

	!### Array to store the emissivity and neutrino spectra ###
	REAL*8, ALLOCATABLE, DIMENSION(:) :: j, num, anum, nue, anue, ecr

	!### Pointer to the emission frequencies ###
	REAL*8, POINTER :: nu(:)

	!### Useful Integers ###
	INTEGER :: nphbins, ncrbins

	!### Flag to see if proton-proton interaction is on ###
	LOGICAL :: pp=.FALSE.

END TYPE

CONTAINS

!### Initializes the proton-proton interaction Object ###
SUBROUTINE PP_Init(ppo,physo,nu,nphbins,elo,ehi,ncrbins)
	TYPE(PP), INTENT(INOUT) :: ppo
	TYPE(Phys), INTENT(IN) :: physo

	REAL*8, TARGET :: nu(nphbins)
	REAL*8, INTENT(IN) :: elo, ehi
	INTEGER, INTENT(IN) :: nphbins, ncrbins

	if (ncrbins .LT. 1) then
		write(*,*) 'PP ERROR: The number of Secondary Bins must be larger than 0'
		STOP
	endif

	if (elo .GE. ehi .AND. ncrbins .GE. 2) then
		write(*,*) 'PP ERROR: The lower energy bound must be less than the upper bound'
		STOP
	endif

	if (ncrbins .EQ. 1 .AND. elo .NE. ehi) then
		write(*,*) 'PP WARNING: For a single energy bin, the lower bound is used'
	endif

	!### Allocate space for the emissivity, and secondaries ###
	ALLOCATE(ppo%j(nphbins),ppo%num(ncrbins),ppo%anum(ncrbins),ppo%nue(ncrbins),ppo%anue(ncrbins),ppo%ecr(ncrbins))

	!### Store number of emission bins ###
	ppo%nphbins=nphbins

	!### Store number of Secondary CR bins ###
	ppo%ncrbins=ncrbins
	ppo%po%nbins=ncrbins
	ppo%eo%nbins=ncrbins

	!### Allocate Space for Secondary Electrons and Positrons ###
	CALL Info_InitCR(ppo%eo)
	CALL Info_InitCR(ppo%po)

	!### Initialize Electrons and Positrons and Secondary energy range. ###
	CALL PP_SecondaryInit(physo,ppo%eo,ppo%po,ppo%ecr,elo,ehi,ncrbins)

	!### Set up pointer to frequency ###
	ppo%nu => nu

END SUBROUTINE PP_Init

!### Initializes CR Secondaries ###
SUBROUTINE PP_SecondaryInit(physo,eo,po,ecr,elo,ehi,nbins)
	TYPE(CR), INTENT(INOUT) :: eo
	TYPE(CR), INTENT(INOUT) :: po
	TYPE(Phys), INTENT(IN) :: physo

	REAL*8, INTENT(INOUT) :: ecr(nbins)
	REAL*8, INTENT(IN) :: elo, ehi
	REAL*8 :: dloge
	INTEGER, INTENT(IN) :: nbins
	INTEGER :: i

	!### First the electrons ###
	!### Initialize type of particle ###
	eo%q=-physo%q	!Charge of Electron, esu
	eo%m=physo%me	!Mass of Electron, g

	!### Initialize the kinetic energy in ergs ###
	dloge=(dlog(ehi)-dlog(elo))/REAL(nbins-1,8)

	!### Now we will set the energies ###
	do i=1,nbins
		eo%e(i)=dlog(elo)+REAL(i-1,8)*dloge
	enddo

	eo%e=dexp(eo%e)

	!### Second the positrons ###
	!### Initialize type of particle ###
	po%q=physo%q	!Charge of Positron, esu
	po%m=physo%me	!Mass of Electron, g

	!### Now the energies
	po%e=eo%e

	!### Finally the other energy range ###
	ecr=eo%e

	!### Calculate useful quantities for the electrons and positrons. ###
	CALL Info_CRCalc(physo,eo)
	CALL Info_CRCalc(physo,po)

END SUBROUTINE PP_SecondaryInit

!### Clean up proton-proton interaction Object ###
SUBROUTINE PP_Destroy(ppo)
	TYPE (PP), INTENT(INOUT) :: ppo

	!### Clean-up pointers ###
	NULLIFY(ppo%nu)

	!### De-Allocate Arrays ###
	DEALLOCATE(ppo%j,ppo%num,ppo%anum,ppo%nue,ppo%anue,ppo%ecr)

	!### De-Allocate Secondaries ###
	CALL Info_DestroyCR(ppo%eo)
	CALL Info_DestroyCr(ppo%po)

END SUBROUTINE PP_Destroy

!### Does the proton-proton interaction calculation ###
!### The formalism in Kelner et.al. PhysRev. D 74 (2006) 034018 is used through out (hence forth refered to as K06). ###
!### Additionally the parameterization defined in Kamae et. al. ApJ 647 (2006) 692 is used for the low energy cross-section (E<100 GeV) ###
!### Note that Kamae (2006) has an errata section that fixes several of the equations. ###
!### Psuedo units means that those are the units if the constants were included.  The constants to be included are found at the end of the problem. ###
SUBROUTINE PP_Calc(ppo,cro,physo,np,nH)
	TYPE(PP), INTENT(INOUT) :: ppo
	TYPE(CR), INTENT(IN) :: cro
	TYPE(Phys), INTENT(IN) :: physo

	REAL*8, INTENT(IN) :: np, nH
	REAL*8 :: npi(cro%nbins), epi(cro%nbins)
	REAL*8 :: jl, ju, eph, emin, eth, Kpi, n, L, slope, jle, jue, jlnue, junue, jlnum, junum, Pp, ntilde, epil, epih
	INTEGER :: i, j
	LOGICAL :: prevpass

	!### Going to integrate it anyways so might as well zero it now. ###
	ppo%j=0d0
	ppo%num=0d0
	ppo%anum=0d0
	ppo%nue=0d0
	ppo%anue=0d0
	ppo%eo%n=0d0
	ppo%po%n=0d0

	!### number of protons and Hydrogen ###
	n=np+nH

	!### We don't need to run this routine if the density is zero or if the CR's are not protons ###
	if (n .EQ. 0d0 .OR. cro%m .NE. physo%mp) then
		RETURN
	endif

	!### We will calculate the photons generated from pion production first ###
	!### Pion Production Threshold Value, 1.22 GeV ###
	eth=1.22d9*physo%eV

	!### Other useful factors ###
	Kpi=0.17d0
	prevpass=.FALSE.
	ntilde=1d0

	!### Using Eq. 77 we will construct a pion spectrum to be used for the low energy end of the gamma ray spectrum ###
	!### This excludes multiplicative factors which will be multiplied back in at the end. ###
	do i=1,cro%nbins
		!### Eq. 11 of K06 ###
		L=dlog(cro%te(i)/(1d12*physo%eV))

		!### Particle momentum in GeV/c ###
		Pp=dsqrt(cro%te(i)**2/physo%c**2-cro%m**2*physo%c**2)/(1d9*physo%eV)*physo%c

		!### Check to see if we are above threshold ###
		if (cro%te(i) .GE. eth) then
			!### Eq. 77  in K06 ###
			epi(i)=Kpi*(cro%te(i)-physo%mpi0*physo%c**2)
			npi(i)=sigmai(L,cro%te(i),eth,Pp,physo%eV)*cro%n(i)*cro%e(i)/cro%te(i)
		else
			!### Otherwise zero out everything ###
			epi(i)=Kpi*(eth-physo%mpi0*physo%c**2)
			npi(i)=0d0
		endif
	enddo

	!### Marches through emission frequencies backwards.  We need to do this in order to match the low energy spectrum with the high energy spectrum ###
	do j=ppo%nphbins,1,-1

		!### Calculate the photon energy and the minimum energy for the integral to be done. ###
		!### Found under Eq. 78 in K06, converted for use with Ep rather than Epi ###
		eph=physo%h*ppo%nu(j)

		!### If we are above 100 GeV in photon energy we will treat the problem as per Eq. 71 in K06 ###
		if (eph .GE. 100d9*physo%eV) then

			!### Flag to let us know that the previous pass was above or at 100 GeV ###
			prevpass=.TRUE.

			!### Eq. 11 of K06 ###
			L=dlog(cro%te(1)/(1d12*physo%eV))

			!### Particle momentum in GeV/c ###
			Pp=dsqrt(cro%te(1)**2/physo%c**2-cro%m**2*physo%c**2)/(1d9*physo%eV)*physo%c

			!### We will calculate the first bin in the integral for later use. ###
			!### We need to make sure that the cosmic ray in question can actually produce any radiation at this energy ###
			if (cro%te(1) .GE. eph) then
				!### Eq. 71 of K06 ###
				!### This leaves out constants which will be multiplied back in at the end ###
				!### Psuedo units of photons/cm^3/s/erg/erg ###
				jl=sigmai(L,cro%te(1),eth,Pp,physo%eV)*Fg(eph*cro%ite(1),L)*cro%ite(1)*cro%n(1)
			else
				jl=0d0
			endif

			!### Integrates CR Spectrum ###
			do i=1,cro%nbins-1
				!### Eq. 11 of K06 ###
				L=dlog(cro%te(i+1)/(1d12*physo%eV))

				!### Particle momentum in GeV/c ###
				Pp=dsqrt(cro%te(i+1)**2/physo%c**2-cro%m**2*physo%c**2)/(1d9*physo%eV)*physo%c

				!### We will calculate the first bin in the integral for later use. ###
				!### We need to make sure that the cosmic ray in question can actually produce any radiation at this energy ###
				if (cro%te(i+1) .GE. eph) then
					!### Eq. 71 of K06 ###
					!### This leaves out constants which will be multiplied back in at the end ###
					!### Psuedo units of photons/cm^3/s/erg/erg ###
					ju=sigmai(L,cro%te(i+1),eth,Pp,physo%eV)*Fg(eph*cro%ite(i+1),L)*cro%ite(i+1)*cro%n(i+1)
				else
					ju=0d0
				endif

				!### We will assume that the function produced by the product of the radiation power ###
				!### with the CR spectrum is a peice-wise powerlaw for the sake of integration. ###

				!### Thus first we will calculate the powerlaw slope ###
				slope=(dlog(MAX(ju,1d-100))-dlog(MAX(jl,1d-100)))/(dlog(cro%e(i+1))-dlog(cro%e(i)))

				!### Now we integrate over cosmic ray energy using the assumed powerlaw ###
				!### We need to test the slope to see if it is close to -1. ###
				!### If it is close to -1 then an edge case occurs and we need to change the integration ###
				!### ppo%j should be in psuedo-units of photons/cm^3/s/erg ###
				!### This is excluding the normalization constants at this point which will put it into proper units. ###
				if (slope .LT. -.999999d0 .AND. slope .GT. -1.000001d0) then
					!### We will now calculate the normalization of the powerlaw ###
					if (jl .NE. 0d0) then
						jl=jl/cro%e(i)**slope
					else
						jl=ju/cro%e(i+1)**slope
					endif
	
					ppo%j(j)=jl*dlog(cro%e(i+1)/cro%e(i))+ppo%j(j)
				else
					ppo%j(j)=1d0/(slope+1d0)*(ju*cro%e(i+1)-jl*cro%e(i))+ppo%j(j)
				endif

				!### Store the upper value for the next update. ###
				jl=ju
			enddo
		else
			!### We are below 100 GeV in photon energy so we will now use Eq. 78 instead to derive the photon spectrum ###

			!### First we will check to see if the previous pass was at or above 100 GeV ###
			if (prevpass) then
				!### Set flag to false ###
				prevpass=.FALSE.
				
				!### We will now calculate the normalization, value of ntilde, needed to match the high energy spectrum to the low energy one. ###
				ntilde=0d0

				!### First we need to set the lower bound for the integral, this is found in the line below Eq. 78 of K06 ###
				emin=physo%h*ppo%nu(j+1)+0.25d0*(physo%mpi0*physo%c**2)**2/(physo%h*ppo%nu(j+1))

				!### Now we will calculate the first value for the integral ###
				if (epi(1) .GE. emin) then
					!### Eq. 78 in K06, excluding numberi ###
					jl=npi(1)/dsqrt(epi(1)**2-(physo%mpi0*physo%c**2)**2)
				else
					jl=0d0
				endif

				!### We will do the rest of the integral now. ###
				do i=1,cro%nbins-1
					if (epi(i+1) .GE. emin) then
						!### Eq. 78 in K06, excluding numberi ###
						ju=npi(i+1)/dsqrt(epi(i+1)**2-(physo%mpi0*physo%c**2)**2)
					else
						ju=0d0
					endif

					!### We will assume that the function produced by the product of the radiative power ###
					!### with the CR spectrum is a peice-wise powerlaw for the sake of integration. ###

					!### Safety check as the energy below threshold has the same value ###
					if (epi(i+1) .NE. epi(i)) then

						epil=epi(i)
						epih=epi(i+1)

						!### This will adjust the lower bound in the case that the lower bound of the integral falls in between zones ###
						if (epi(i) .LT. emin .AND. epi(i+1) .GE. emin) then
							epil=emin
							slope=(dlog(MAX(npi(i+1),1d-100))-dlog(MAX(npi(i),1d-100)))/(dlog(epi(i+1))-dlog(epi(i)))
							jl=npi(i)*(emin/epi(i))**slope/dsqrt(emin**2-(physo%mpi0*physo%c**2)**2)
						endif

						!### Thus first we will calculate the powerlaw slope ###
						slope=(dlog(MAX(ju,1d-100))-dlog(MAX(jl,1d-100)))/(dlog(epih)-dlog(epil))

						!### Now we integrate over cosmic ray energy using the assumed powerlaw ###
						!### We need to test the slope to see if it is close to -1. ###
						!### If it is close to -1 then an edge case occurs and we need to change the integration ###
						!### ppo%j should be in psuedo-units of photons/cm^3/s/erg ###
						!### This is excluding the normalization constants at this point which will put it into proper units. ###
						if (slope .LT. -.999999d0 .AND. slope .GT. -1.000001d0) then
							!### We will now calculate the normalization of the powerlaw ###
							if (jl .NE. 0d0) then
								jl=jl/epil**slope
							else
								jl=ju/epih**slope
							endif

							ntilde=jl*dlog(epih/epil)+ntilde
						else
							ntilde=1d0/(slope+1d0)*(ju*epih-jl*epil)+ntilde
						endif

					endif

					!### Store the upper value for the next update. ###
					jl=ju
				enddo

				!### Now we will complete our normalization by inverting Eq. 77 and 78 and solving for ntilde ###
				if (ntilde .EQ. 0d0) then
					ntilde=0d0
				else
					ntilde=ppo%j(j+1)*Kpi/(2d0*ntilde)
				endif
			endif

			!### Now that we have completed the normalization we can move on to actually computing Eq. 78 properly ###
			!### First we need to set the lower bound for the integral, this is found in the line below Eq. 78 of K06 ###
			emin=physo%h*ppo%nu(j)+0.25d0*(physo%mpi0*physo%c**2)**2/(physo%h*ppo%nu(j))

			!### Now we will calculate the first value for the integral ###
			if (epi(1) .GE. emin) then
				!### Eq. 78 in K06, excluding numberi ###
				jl=npi(1)/dsqrt(epi(1)**2-(physo%mpi0*physo%c**2)**2)
			else
				jl=0d0
			endif

			!### We will do the rest of the integral now. ###
			do i=1,cro%nbins-1
				if (epi(i+1) .GE. emin) then
					!### Eq. 78 in K06, excluding numberi ###
					ju=npi(i+1)/dsqrt(epi(i+1)**2-(physo%mpi0*physo%c**2)**2)
				else
					ju=0d0
				endif

				!### We will assume that the function produced by the product of the radiative power ###
				!### with the CR spectrum is a peice-wise powerlaw for the sake of integration. ###

				!### Safety check as the energy below threshold has the same value ###
				if (epi(i+1) .NE. epi(i)) then

					epil=epi(i)
					epih=epi(i+1)

					!### This will adjust the lower bound in the case that the lower bound of the integral falls in between zones ###
					if (epi(i) .LT. emin .AND. epi(i+1) .GE. emin) then
						epil=emin
						slope=(dlog(MAX(npi(i+1),1d-100))-dlog(MAX(npi(i),1d-100)))/(dlog(epi(i+1))-dlog(epi(i)))
						jl=npi(i)*(emin/epi(i))**slope/dsqrt(emin**2-(physo%mpi0*physo%c**2)**2)
					endif

					!### Thus first we will calculate the powerlaw slope ###
					slope=(dlog(MAX(ju,1d-100))-dlog(MAX(jl,1d-100)))/(dlog(epih)-dlog(epil))

					!### Now we integrate over cosmic ray energy using the assumed powerlaw ###
					!### We need to test the slope to see if it is close to -1. ###
					!### If it is close to -1 then an edge case occurs and we need to change the integration ###
					!### ppo%j should be in psuedo-units of photons/cm^3/s/erg ###
					!### This is excluding the normalization constants at this point which will put it into proper units. ###
					if (slope .LT. -.999999d0 .AND. slope .GT. -1.000001d0) then
						!### We will now calculate the normalization of the powerlaw ###
						if (jl .NE. 0d0) then
							jl=jl/epil**slope
						else
							jl=ju/epih**slope
						endif

						ppo%j(j)=jl*dlog(epih/epil)+ppo%j(j)
					else
						ppo%j(j)=1d0/(slope+1d0)*(ju*epih-jl*epil)+ppo%j(j)
					endif
				endif

				!### Store the upper value for the next update. ###
				jl=ju
			enddo

			!### Now we will complete Eq. 77 and 78 in K06 ###
			ppo%j(j)=2d0*ntilde/Kpi*ppo%j(j)
		endif
	enddo

	!### Now we will multiply in the constants we left out of Eq. 71 and Eq. 78. ###
	!### Plus we will convert the distribution to units of emissivity which are ergs/s/cm^3/str/Hz ###
	ppo%j=n*physo%c*physo%h**2*ppo%nu*0.25d0/physo%pi*ppo%j

	!### Handles the Secondary CR production. ###

	!### Marches through emission frequencies ###
	do j=1,ppo%ncrbins

		!### Eq. 11 of K06 ###
		L=dlog(cro%te(1)/(1d12*physo%eV))

		!### Particle momentum in GeV/c ###
		Pp=dsqrt(cro%te(1)**2/physo%c**2-cro%m**2*physo%c**2)/(1d9*physo%eV)*physo%c

		!### Calculate the photon source function for the first bin for use later ###
		!### Electrons first ###
		if (cro%te(1) .LT. ppo%eo%te(j)) then
			!### If below the minimum energy then there is no way to produce electrons of this energy ###
			jle=0d0
		else
			!### Eq. 71 of K06 ###
			!### This leaves out constants which will be multiplied back in at the end ###
			!### Psuedo units of particles/cm^3/s/erg/erg ###
			!### Appears to be typo in K06, x=E_e/E_p not E_pi ###
			jle=sigmai(L,cro%te(1),eth,Pp,physo%eV)*Fe(ppo%eo%te(j)*cro%ite(1),L)*cro%ite(1)*cro%n(1)

		endif

		!### Now neutrinos ###
		if (cro%te(1) .LT. ppo%ecr(j)) then
			!### If below the minimum energy then there is no way to produce neutrinos of this energy ###
			jlnue=0d0
			jlnum=0d0
		else
			!### Eq. 71 of K06 ###
			!### This leaves out constants which will be multiplied back in at the end ###
			!### Psuedo units of particles/cm^3/s/erg/erg ###
			!### electron neutrinos ###
			jlnue=sigmai(L,cro%te(1),eth,Pp,physo%eV)*Fe(ppo%ecr(j)*cro%ite(1),L)*cro%ite(1)*cro%n(1)

			!### muon neutrinos ###
			jlnum=sigmai(L,cro%te(1),eth,Pp,physo%eV)*(Fe(ppo%ecr(j)*cro%ite(1),L)+Fnum(ppo%ecr(j)*cro%ite(1),L))*cro%ite(1)*cro%n(1)
		endif

		!### Integrates CR Spectrum ###
		do i=1,cro%nbins-1
			!### Eq. 11 of K06 ###
			L=dlog(cro%te(i+1)/(1d12*physo%eV))

			!### Particle momentum in GeV/c ###
			Pp=dsqrt(cro%te(i+1)**2/physo%c**2-cro%m**2*physo%c**2)/(1d9*physo%eV)*physo%c

			!### Calculate the photon source function ###
			!### Electrons first ###
			if (cro%te(i+1) .LT. ppo%eo%te(j)) then
				!### If below the minimum energy then there is no way to produce electrons of this energy ###
				jue=0d0
			else
				!### Eq. 71 of K06 ###
				!### This leaves out constants which will be multiplied back in at the end ###
				!### Psuedo units of particles/cm^3/s/erg/erg ###
				!### Appears to be typo in K06, x=E_e/E_p not E_pi ###
				jue=sigmai(L,cro%te(i+1),eth,Pp,physo%eV)*Fe(ppo%eo%te(j)*cro%ite(i+1),L)*cro%ite(i+1)*cro%n(i+1)

			endif

			!### Now neutrinos ###
			if (cro%te(i+1) .LT. ppo%ecr(j)) then
				!### If below the minimum energy then there is no way to produce neutrinos of this energy ###
				junue=0d0
				junum=0d0
			else
				!### Eq. 71 of K06 ###
				!### This leaves out constants which will be multiplied back in at the end ###
				!### Psuedo units of particles/cm^3/s/erg/erg ###
				!### electron neutrinos ###
				junue=sigmai(L,cro%te(i+1),eth,Pp,physo%eV)*Fe(ppo%ecr(j)*cro%ite(i+1),L)*cro%ite(i+1)*cro%n(i+1)

				!### muon neutrinos ###
				junum=sigmai(L,cro%te(i+1),eth,Pp,physo%eV)*(Fe(ppo%ecr(j)*cro%ite(i+1),L)+Fnum(ppo%ecr(j)*cro%ite(i+1),L))*cro%ite(i+1)*cro%n(i+1)
			endif

			!### We will assume that the function produced by the product of the synchrotron power ###
			!### with the CR spectrum is a peice-wise powerlaw for the sake of integration. ###

			!### Electrons first ###
			!### Thus first we will calculate the powerlaw slope ###
			slope=(dlog(MAX(jue,1d-100))-dlog(MAX(jle,1d-100)))/(dlog(cro%e(i+1))-dlog(cro%e(i)))

			!### Now we integrate over cosmic ray energy using the assumed powerlaw ###
			!### We need to test the slope to see if it is close to -1. ###
			!### If it is close to -1 then an edge case occurs and we need to change the integration ###
			!### ppo%j should be in psuedo-units of photons/cm^3/s/erg ###
			!### This is excluding the normalization constants at this point which will put it into proper units. ###
			if (slope .LT. -.999999d0 .AND. slope .GT. -1.000001d0) then
				!### We will now calculate the normalization of the powerlaw ###
				if (jle .NE. 0d0) then
					jle=jle/cro%e(i)**slope
				else
					jle=jue/cro%e(i+1)**slope
				endif

				ppo%eo%n(j)=jle*dlog(cro%e(i+1)/cro%e(i))+ppo%eo%n(j)
			else
				ppo%eo%n(j)=1d0/(slope+1d0)*(jue*cro%e(i+1)-jle*cro%e(i))+ppo%eo%n(j)
			endif

			!### Store the upper value for the next update. ###
			jle=jue

			!### Electron neutrinos next ###
			!### Thus first we will calculate the powerlaw slope ###
			slope=(dlog(MAX(junue,1d-100))-dlog(MAX(jlnue,1d-100)))/(dlog(cro%e(i+1))-dlog(cro%e(i)))

			!### Now we integrate over cosmic ray energy using the assumed powerlaw ###
			!### We need to test the slope to see if it is close to -1. ###
			!### If it is close to -1 then an edge case occurs and we need to change the integration ###
			!### ppo%j should be in psuedo-units of photons/cm^3/s/erg ###
			!### This is excluding the normalization constants at this point which will put it into proper units. ###
			if (slope .LT. -.999999d0 .AND. slope .GT. -1.000001d0) then
				!### We will now calculate the normalization of the powerlaw ###
				if (jle .NE. 0d0) then
					jlnue=jlnue/cro%e(i)**slope
				else
					jlnue=junue/cro%e(i+1)**slope
				endif

				ppo%nue(j)=jlnue*dlog(cro%e(i+1)/cro%e(i))+ppo%nue(j)
			else
				ppo%nue(j)=1d0/(slope+1d0)*(junue*cro%e(i+1)-jlnue*cro%e(i))+ppo%nue(j)
			endif

			!### Store the upper value for the next update. ###
			jlnue=junue

			!### Muon neutrinos last ###
			!### Thus first we will calculate the powerlaw slope ###
			slope=(dlog(MAX(junum,1d-100))-dlog(MAX(jlnum,1d-100)))/(dlog(cro%e(i+1))-dlog(cro%e(i)))

			!### Now we integrate over cosmic ray energy using the assumed powerlaw ###
			!### We need to test the slope to see if it is close to -1. ###
			!### If it is close to -1 then an edge case occurs and we need to change the integration ###
			!### ppo%j should be in psuedo-units of photons/cm^3/s/erg ###
			!### This is excluding the normalization constants at this point which will put it into proper units. ###
			if (slope .LT. -.999999d0 .AND. slope .GT. -1.000001d0) then
				!### We will now calculate the normalization of the powerlaw ###
				if (jle .NE. 0d0) then
					jlnum=jlnum/cro%e(i)**slope
				else
					jlnum=junum/cro%e(i+1)**slope
				endif

				ppo%num(j)=jlnum*dlog(cro%e(i+1)/cro%e(i))+ppo%num(j)
			else
				ppo%num(j)=1d0/(slope+1d0)*(junum*cro%e(i+1)-jlnum*cro%e(i))+ppo%num(j)
			endif

			!### Store the upper value for the next update. ###
			jlnum=junum
		enddo
	enddo

	!### Now we will multiply in the constants we left out of Eq. 71. ###
	!### Units will be particles/cm^3/s/erg ###
	ppo%eo%n=n*physo%c*ppo%eo%n
	ppo%nue=n*physo%c*ppo%nue
	ppo%num=n*physo%c*ppo%num

	!### The anti-particles are symmetric ###
	ppo%po%n=ppo%eo%n
	ppo%anue=ppo%nue
	ppo%anum=ppo%num

END SUBROUTINE PP_Calc

!### Inclusive Cross-section in cm^2, Eq. 79 in K06 and Kamae (2006) ###
REAL*8 FUNCTION sigmai(L,Ep,Eth,Pp,eV)
	REAL*8 :: L, Ep, Eth, Pp, eV, x, Epgev
	REAL*8 :: snd, sdiff, sdel, sres
	REAL*8 :: a0, a1, a2, a3, a4, a5, a6, a7
	REAL*8 :: b0, b1, c0, c1, c2
	REAL*8 :: d0, d1, d2, d3, d4, d5, d6
	REAL*8 :: e0, e1, f0, f1, f2, f3, f4
	REAL*8 :: g0, g1, g2, g3, g4


	!### If the momentum is lower than 100 GeV/c use the cross-section defined in Kamae (2006), otherwise use Eq. 79 in K06 ###
	if (Pp .LT. 1d2) then
		!### Log of the momentum ###
		x=dlog10(Pp)

		!### Proton energy in GeV ###
		Epgev=Ep/(1d9*eV)

		!### Constants defined in Table 1 of Kamae (2006) Errata ###
		a0=0.1176d0
		a1=0.3829d0
		a2=23.10d0
		a3=6.454d0
		a4=-5.764d0
		a5=-23.63d0
		a6=94.75d0
		a7=0.02667d0

		b0=11.34d0
		b1=23.72d0

		c0=28.5d0
		c1=-6.133d0
		c2=1.464d0

		!### This is a bug fix.  Apparently this nees to be exactly log10 of 2.25 rather than what Kamae has of 0.3522 ###
		d0=dlog10(2.25d0)
		d1=0.1530d0
		d2=1.498d0
		d3=2.0d0
		d4=30.0d0
		d5=3.155d0
		d6=1.042d0

		e0=5.922d0
		e1=1.632d0

		f0=0.0834d0
		f1=9.5d0
		f2=-5.5d0
		f3=1.68d0
		f4=3134d0

		g0=0.0004257d0
		g1=4.5d0
		g2=-7.0d0
		g3=2.1d0
		g4=503.5d0


		!### Non-diffractive cross-section, Eq. 1 in Kamae (2006) Errata ###
		if (Pp .LT. 1d0) then
			snd=0d0
		elseif (Pp .GE. 1d0 .AND. Pp .LT. 1.3d0) then
			snd=0.57d0*(x/a0)**(1.2d0)*(a2+a3*x**2+a4*x**3+a5*dexp(-a6*(x+a7)**2))
		elseif (Pp .GE. 1.3d0 .AND. Pp .LT. 2.4d0) then
			snd=(b0*abs(a1-x)+b1*abs(a0-x))/(a1-a0)
		elseif (Pp .GE. 2.4d0 .AND. Pp .LT. 10d0) then
			snd=a2+a3*x**2+a4*x**3+a5*dexp(-a6*(x+a7)**2)
		else
			snd=c0+c1*x+c2*x**2
		endif

		!### Diffractive cross-section, Eq. 2 in Kamae (2006) ###
		if (Pp .LT. 2.25d0) then
			sdiff=0d0
		elseif (Pp .GE. 2.25d0 .AND. Pp .LT. 3.2d0) then
			sdiff=dsqrt((x-d0)/d1)*(d2+d3*dlog10(d4*(x-0.25d0))+d5*x**2-d6*x**3)
		elseif (Pp .GE. 3.2d0 .AND. Pp .LT. 100d0) then
			sdiff=d2+d3*dlog10(d4*(x-0.25d0))+d5*x**2-d6*x**3
		else
			sdiff=e0+e1*x
		endif


		!### Delta 1232 resonance cross-section, Eq. 3 in Kamae (2006) Errata ###
		if (Epgev .LT. 1.4d0) then
			sdel=0d0
		elseif (Epgev .GE. 1.4d0 .AND. Epgev .LT. 1.6d0) then
			sdel=f0*Epgev**10
		elseif (Epgev .GE. 1.6d0 .AND. Epgev .LT. 1.8d0) then
			sdel=f1*dexp(-f2*(Epgev-f3)**2)
		elseif (Epgev .GE. 1.8d0 .AND. Epgev .LT. 10d0) then
			sdel=f4*Epgev**(-10)
		else
			sdel=0d0
		endif

		!### 1600 Resonance cross-section, Eq. 4 in Kamae (2006) Errata ###
		if (Epgev .LT. 1.6d0) then
			sres=0d0
		elseif (Epgev .GE. 1.6d0 .AND. Epgev .LT. 1.9d0) then
			sres=g0*Epgev**14
		elseif (Epgev .GE. 1.9d0 .AND. Epgev .LT. 2.3d0) then
			sres=g1*dexp(-g2*(Epgev-g3)**2)
		elseif (Epgev .GE. 2.3d0 .AND. Epgev .LT. 20d0) then
			sres=g4*Epgev**(-6)
		else
			sres=0d0
		endif

		sigmai=(snd+sdiff+sdel+sres)*1d-27
	else
		!### Eq. 79 in K06 ###
		sigmai=(34.3d0+1.88d0*L+0.25d0*L**2)*(1d0-(Eth/Ep)**4)**2*1d-27
	endif

END FUNCTION sigmai

!### Eq. 58 of K06 ###
REAL*8 FUNCTION Fg(x,L)
	REAL*8 :: x, L, B, beta, k

	!### Eq. 59, 60, 61 of K06 ###
	B=1.3d0+0.14d0*L+0.011d0*L**2
	beta=1d0/(1.79d0+0.11d0*L+0.008d0*L**2)
	k=1d0/(0.801d0+0.049d0*L+0.014*L**2)

	!### Eq. 58 of K06
	Fg=B*dlog(x)/x*((1d0-x**beta)/(1d0+k*x**beta*(1d0-x**beta)))**4 &
	*(1d0/dlog(x)-4d0*beta*x**beta/(1d0-x**beta)-4d0*k*beta*x**beta*(1d0-2d0*x**beta)/(1d0+k*x**beta*(1d0-x**beta)))

END FUNCTION Fg

!### Eq. 62 of K06 ###
REAL*8 FUNCTION Fe(x,L)
	REAL*8 :: x, L, B, beta, k

	!### Eq. 63, 64, 65 of K06 ###
	B=1d0/(69.5d0+2.65d0*L+0.3d0*L**2)
	beta=1d0/(0.201d0+0.062d0*L+0.00042d0*L**2)**(0.25d0)
	k=(0.279d0+0.141*L+0.0172*L**2)/(0.3d0+(2.3d0+L)**2)

	!### Eq. 62 of K06 ###
	Fe=B*(1d0+k*dlog(x)**2)**3/(x*(1d0+0.3d0/x**beta))*(-dlog(x))**5

END FUNCTION Fe

!### Eq. 66 of K06 ###
REAL*8 FUNCTION Fnum(x,L)
	REAL*8 :: x, L, B, beta, k, y

	y=x/0.427d0

	!### Spectrum cuts off sharply at this point. ###
	if (y .GT. 1d0) then
		Fnum=0d0
		RETURN
	endif

	!### Eq. 67, 68, 69 of K06 ###
	B=1.75d0+0.204d0*L+0.01d0*L**2
	beta=1d0/(1.67d0+0.111d0*L+0.0038d0*L**2)
	k=1.07d0-0.086d0*L+0.002*L**2

	!### Eq. 66 of K06
	Fnum=B*dlog(y)/y*((1d0-y**beta)/(1d0+k*y**beta*(1d0-y**beta)))**4 &
	*(1d0/dlog(y)-4d0*beta*y**beta/(1d0-y**beta)-4d0*k*beta*y**beta*(1d0-2d0*y**beta)/(1d0+k*y**beta*(1d0-y**beta)))

END FUNCTION Fnum

END MODULE Mod_PP
