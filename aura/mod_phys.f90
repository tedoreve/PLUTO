MODULE Mod_Phys
!############################################################################
!#  This module contains the definition of the physical units stored in
!#  the phys object.
!############################################################################

IMPLICIT NONE

TYPE :: Phys
	REAL*8 c, q, me, r0, pi, sigmat, eV, h, kb, alpha, mp, mn, mpi0, mpip, mpin, mmu
END TYPE

CONTAINS

!### This subroutine initializes the physical units ###
SUBROUTINE Phys_Init(physo)
	TYPE(Phys) :: physo

	!### Physical Constants in cgs ###
	physo%c=2.99792458d10				!Speed of Light, cm/s
	physo%q=4.8032068d-10 				!Charge, esu
	physo%me=9.1093897d-28 				!Mass of electron, g
	physo%r0=physo%q**2/(physo%me*physo%c**2) 	!Classical Electron Radius, cm
	physo%pi=4.0d0*datan(1.0d0) 			!pi
	physo%sigmat=8d0*physo%pi/3d0*physo%r0**2 	!Thomson Cross Section cm^2
	physo%eV=1.60217646d-12 			!Number of ergs in an electron-Volt (eV)
	physo%h=6.6260755d-27 				!Planck Constant, erg*s
	physo%kb=1.380658d-16				!Boltzmann constant, erg/K
	physo%alpha=7.29735308d-3 			!Fine Structure Constant
	physo%mp=1.6726231d-24 				!Mass of proton, g
	physo%mn=1.6749286d-24 				!Mass of neutron, g
	physo%mpi0=2.406119d-25 			!Mass of neutral pion, g
	physo%mpip=2.488017d-25 			!Mass of pi^+, g
	physo%mpin=physo%mpip 				!Mass of pi^-, g
	physo%mmu=1.8835317d-25 			!Mass of muon, g

END SUBROUTINE Phys_Init

END MODULE Mod_Phys
