MODULE Mod_User

!############################################################################
!#  This section contains the user defined CR spectrum and environoment
!############################################################################

!### Modules to be Included ###
USE Mod_Info

IMPLICIT NONE

CONTAINS

!### Initializes User Defined Parameters ###
SUBROUTINE User_Init(info)
	TYPE(Infoc) :: info

	if (info%verbose) write(*,*) 'Defining User Parameters'

	!### Sets up environment ###
	CALL User_Env(info%envo,info%physo)

	!### Sets up particle spectrum ###
	CALL User_CR(info%cro,info%physo)

END SUBROUTINE User_Init

!### Initializes environment in cgs units.  Photon spectra should be stored in units of [photons/cm^3/erg] ###
SUBROUTINE User_Env(envo,physo)
	TYPE(Env) :: envo
	TYPE(Phys) :: physo

	REAL*8 :: dlognu, nulo, nuhi, Tbb
	INTEGER :: i

	envo%ne=1d0			!Electron Number Density, cm^{-3}
	envo%np=1d0			!Proton Number Density, cm^{-3}
	envo%nalpha=0d0			!Alpha Particle Number Density, cm^{-3}, Bremsstrahlung only
	envo%nH=0d0			!Neutral Hydrogen Number Density, cm^{-3}
	envo%nHe=0d0			!Neutral Helium Number Density, cm^{-3}, Bremsstrahlung only
	envo%B=1d-6			!Magnetic Field Strength, G
	envo%theta=0.5d0*physo%pi 	!Angle of magnetic field with respect to the observer, radians, for SED calculations this should be set to pi/2

	!### Setup ambient photon field. ###
	!### First we will set up the frequency range in Hz. ###
	nulo=1d0
	nuhi=1d14

	!### Set size of logarithmic spacing ###
	dlognu=(dlog(nuhi)-dlog(nulo))/REAL(envo%nphbins-1,8)

	!### Calculate frequencies in log space ###
	do i=1,envo%nphbins
		envo%nuph(i)=dlog(nulo)+REAL(i-1,8)*dlognu
	enddo

	!### Convert back to normal space ###
	envo%nuph=dexp(envo%nuph)

	!### Now we will calculate the ambient photon field ###
	!### This should be in units of photons/cm^3/erg ###
	!### Current set up for single temperature blackbody but it can handle any arbitrary isotropic photon field ###
	Tbb=2.73d0	!Blackbody Temperature, K

	do i=1,envo%nphbins
		envo%nph(i)=(physo%pi**2*(0.5d0*physo%h/physo%pi*physo%c)**3)**(-1) &
		*(physo%h*envo%nuph(i))**2/(dexp(physo%h*envo%nuph(i)/(physo%kb*Tbb))-1d0)
	enddo

END SUBROUTINE User_Env

!### Initializes particle spectrum and particle type in cgs units.  Particle spectra should be stored in units of [particles/cm^3/erg] and energy should be in ergs ###
!### Note that the particle spectrum is dN/dE where E is the kinetic energy, same with the stored CR energy.  DO NOT USE TOTAL ENERGY ### 
SUBROUTINE User_CR(cro,physo)
	TYPE(CR) :: cro
	TYPE(Phys) :: physo

	REAL*8 :: dloge, emin, emax
	INTEGER :: i

	!### Initialize type of particle ###
	cro%q=-physo%q	!Charge of CR particle, esu
	cro%m=physo%me	!Mass of CR particle, g

	!### Initialize the kinetic energy in ergs ###
	!### First we will set the logarithmic spacing ###
	emin=1d3*physo%eV
	emax=1d15*physo%eV

	dloge=(dlog(emax)-dlog(emin))/REAL(cro%nbins-1,8)

	!### Now we will set the energies ###
	do i=1,cro%nbins
		cro%e(i)=dlog(emin)+REAL(i-1,8)*dloge
	enddo

	cro%e=dexp(cro%e)

	!### Now we will set the spectrum in particles/cm^3/erg ###
	do i=1,cro%nbins
		cro%n(i)=(cro%e(i)/(cro%m*physo%c**2)+1d0)**(-2)
	enddo

END SUBROUTINE User_CR

END MODULE Mod_User
