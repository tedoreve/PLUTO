MODULE Mod_IO

!############################################################################
!#  This section contains the IO routines for Aura
!############################################################################

!### Modules to be Included ###
USE Mod_Emis
USE Mod_Info

IMPLICIT NONE

CONTAINS

!### This routine will output the results to text files. ###
!### Current output units are j=[ergs/cm^3/s/Hz/str], n=[particles/cm^3/s/erg] ###
SUBROUTINE IO_Emis(emiso,info)
	TYPE(Emis), INTENT(IN) :: emiso
	TYPE(Infoc), INTENT(IN) :: info

	INTEGER i
	CHARACTER*20 ofile


	!### Synchrotron ###
	if (emiso%synco%synchrotron) then
		if (info%verbose) write(*,*) 'Writing Synchrotron Data'

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A8)') info%oroot, 'sync.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")


		do i=1,emiso%nphbins
			write(12,*) log10(emiso%nuph(i)), log10(MAX(emiso%synco%j(i),1d-100))
		enddo

		CLOSE(12)
	endif

	!### IC ###
	if (emiso%ico%ic) then
		if (info%verbose) write(*,*) 'Writing IC Data'

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A6)') info%oroot, 'ic.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		do i=1,emiso%nphbins
			write(12,*) log10(emiso%nuph(i)), log10(MAX(emiso%ico%j(i),1d-100))
		enddo

		CLOSE(12)
	endif

	!### Bremsstrahlung ###
	if (emiso%bremo%bremsstrahlung) then
		if (info%verbose) write(*,*) 'Writing Bremsstrahlung Data'

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A8)') info%oroot, 'brem.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		do i=1,emiso%nphbins
			write(12,*) log10(emiso%nuph(i)), log10(MAX(emiso%bremo%j(i),1d-100))
		enddo

		CLOSE(12)
	endif

	!### Proton-Proton Interactions ###
	if (emiso%ppo%pp) then
		if (info%verbose) write(*,*) 'Writing Proton-Proton Interaction Data'

		!### Gamma Rays ###

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A7)') info%oroot, 'ppg.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		do i=1,emiso%nphbins
			write(12,*) log10(emiso%nuph(i)), log10(MAX(emiso%ppo%j(i),1d-100))
		enddo

		CLOSE(12)

		!### Electrons ###

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A7)') info%oroot, 'ppe.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		do i=1,emiso%ppo%ncrbins
			write(12,*) log10(emiso%ppo%eo%e(i)), log10(MAX(emiso%ppo%eo%n(i),1d-100))
		enddo

		CLOSE(12)

		!### Electron Neutrinos ###

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A9)') info%oroot, 'ppnue.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		do i=1,emiso%ppo%ncrbins
			write(12,*) log10(emiso%ppo%ecr(i)), log10(MAX(emiso%ppo%nue(i),1d-100))
		enddo

		CLOSE(12)

		!### Muon Neutrinos ###

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A9)') info%oroot, 'ppnum.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		do i=1,emiso%ppo%ncrbins
			write(12,*) log10(emiso%ppo%ecr(i)), log10(MAX(emiso%ppo%num(i),1d-100))
		enddo

		CLOSE(12)
	endif

END SUBROUTINE IO_Emis

!### Testing Routines ###

!### This routine contains useful formulas for testing and will output the results to text files. ###
SUBROUTINE IO_EmisTest(emiso,info)
	TYPE(Emis), INTENT(IN) :: emiso
	TYPE(Infoc), INTENT(IN) :: info

	REAL*8 k, p, a, j, T, nion, sigma, delta
	INTEGER i
	CHARACTER*20 ofile


	!### Synchrotron ###
	if (emiso%synco%synchrotron) then
		if (info%verbose) write(*,*) 'Writing Synchrotron Data'

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A8)') info%oroot, 'sync.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		!### For testing uses Eq. 4.59 from Blumenthal & Gould RMP 42 (1970) 237 ###
!		p=2d0
!		a=0.103d0
!		k=4d0*info%physo%pi*info%physo%me*info%physo%c**2/(8d0*info%physo%pi/3d0)

		do i=1,emiso%nphbins
		!### For testing uses Eq. 4.59 from Blumenthal & Gould RMP 42 (1970) 237 ###
!			j=4d0*info%physo%pi*k*info%physo%q**3*info%envo%B**((p+1d0)/2d0)/(info%physo%me*info%physo%c**2) &
!			*(3d0*info%physo%q/(4d0*info%physo%pi*info%physo%me*info%physo%c))**((p-1d0)/2d0)*a*emiso%nuph(i)**(-(p-1d0)/2d0)

			write(12,*) log10(emiso%nuph(i)), log10(MAX(emiso%synco%j(i),1d-100))!, log10(j)
		enddo

		CLOSE(12)
	endif

	!### IC ###
	if (emiso%ico%ic) then
		if (info%verbose) write(*,*) 'Writing IC Data'

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A6)') info%oroot, 'ic.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		!### For testing uses Eq. 2.65 from Blumenthal & Gould RMP 42 (1970) 237 ###
!		p=2d0
!		a=5.25d0
!		k=info%physo%me*info%physo%c**2/(4d0*info%physo%pi)*info%physo%h
!		T=2.73d0

		do i=1,emiso%nphbins
			!### For testing uses Eq. 2.65 from Blumenthal & Gould RMP 42 (1970) 237 ###
!			j=1d0/info%physo%pi*(info%physo%r0**2*(2d0*info%physo%pi)**3/(info%physo%h**3*info%physo%c**2)) &
!			*k*(info%physo%kb*T)**((p+5d0)/2d0)*a*(info%physo%h*emiso%nuph(i))**(-(p-1d0)/2d0)

			write(12,*) log10(emiso%nuph(i)), log10(MAX(emiso%ico%j(i),1d-100))!, log10(j)
		enddo

		CLOSE(12)
	endif

	!### Bremsstrahlung ###
	if (emiso%bremo%bremsstrahlung) then
		if (info%verbose) write(*,*) 'Writing Bremsstrahlung Data'

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A8)') info%oroot, 'brem.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		!### For testing uses Eq. 3.59 in the weak shielding limit from Blumenthal & Gould RMP 42 (1970) 237 ###
!		p=2d0
!		a=info%physo%h/(4d0*info%physo%pi)!info%physo%me*info%physo%c**2/(4d0*info%physo%pi)
!		nion=info%envo%ne+info%envo%np

		do i=1,emiso%nphbins

			!### For testing uses Eq. 3.59 in the weak shielding limit from Blumenthal & Gould RMP 42 (1970) 237 ###
			!### Weak Shielding (fully ionized case) ###
!			k=info%physo%h*emiso%nuph(i)/(info%physo%me*info%physo%c**2)
!
!			if (k .LT. info%cro%g(1)) then
!				j=info%physo%alpha*info%physo%r0**2*info%physo%c*a/k*nion &
!				*(4d0/3d0*info%cro%g(1)**(-(p-1d0))/(p-1d0)-4d0/3d0*k*info%cro%g(1)**(-p)/p+k**2*info%cro%g(1)**(-(p+1d0))/(p+1d0)) &
!				*4d0*(dlog(2d0*info%cro%g(1)*(info%cro%g(1)-k)/k)-0.5d0)*info%physo%h*emiso%nuph(i)
!			else
!				j=info%physo%alpha*info%physo%r0**2*info%physo%c*a/k*nion &
!				*(4d0/3d0*info%cro%g(1)**(-(p-1d0))/(p-1d0)-4d0/3d0*k*info%cro%g(1)**(-p)/p+k**2*info%cro%g(1)**(-(p+1d0))/(p+1d0)) &
!				*4d0*(dlog(4d0*k)-0.5d0)*info%physo%h*emiso%nuph(i)
!			endif

			write(12,*) log10(emiso%nuph(i)), log10(MAX(emiso%bremo%j(i),1d-100))!, log10(j)
		enddo

		CLOSE(12)
	endif

	!### Proton-Proton Interactions ###
	if (emiso%ppo%pp) then
		if (info%verbose) write(*,*) 'Writing Proton-Proton Interaction Data'

		!### Gamma Rays ###

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A7)') info%oroot, 'ppg.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		!### For Testing using Eq. 23 from Edmon et.al. MNRAS 414 (2011) 3521 ###
!		p=2d0
!		a=info%physo%h/(4d0*info%physo%pi)*(info%physo%mp*info%physo%c**2)**(p-1d0)
!		nion=info%envo%np+info%envo%nH
!		delta=0.14d0*p**(-1.6d0)+0.44d0
!		sigma=3.2d-26*(0.96d0+dexp(4.4d0-2.4*p))

		do i=1,emiso%nphbins
			!### For Testing using Eq. 23 from Edmon et.al. MNRAS 414 (2011) 3521 ###
!			j=a*nion*sigma*info%physo%c/(info%physo%mp*info%physo%c**2)*2d0**(2d0-p)*4d0/(3d0*p)*(info%physo%mpi0/info%physo%mp)**(-p) &
!			*((2d0*info%physo%h*emiso%nuph(i)/(info%physo%mpi0*info%physo%c**2))**delta &
!			+(2d0*info%physo%h*emiso%nuph(i)/(info%physo%mpi0*info%physo%c**2))**(-delta))**(-p/delta)*info%physo%h*emiso%nuph(i)
			write(12,*) log10(emiso%nuph(i)), log10(MAX(emiso%ppo%j(i),1d-100))!, log10(MAX(j,1d-100))
		enddo

		CLOSE(12)

		!### Electrons ###

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A7)') info%oroot, 'ppe.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		do i=1,emiso%ppo%ncrbins
			write(12,*) log10(emiso%ppo%eo%e(i)), log10(MAX(emiso%ppo%eo%n(i),1d-100))
		enddo

		CLOSE(12)

		!### Electron Neutrinos ###

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A9)') info%oroot, 'ppnue.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		do i=1,emiso%ppo%ncrbins
			write(12,*) log10(emiso%ppo%ecr(i)), log10(MAX(emiso%ppo%nue(i),1d-100))
		enddo

		CLOSE(12)

		!### Muon Neutrinos ###

        	!### define the file we'll write to ###
        	WRITE(ofile, '(A5, "_", A9)') info%oroot, 'ppnum.dat'

		!### Open and write to file ###
		OPEN(UNIT=12, FILE=ofile, DEFAULTFILE=info%odump_path, FORM="FORMATTED", STATUS="UNKNOWN")

		do i=1,emiso%ppo%ncrbins
			write(12,*) log10(emiso%ppo%ecr(i)), log10(MAX(emiso%ppo%num(i),1d-100))
		enddo

		CLOSE(12)
	endif

END SUBROUTINE IO_EmisTest

END MODULE Mod_IO
