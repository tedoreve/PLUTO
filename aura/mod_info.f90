MODULE Mod_Info

!############################################################################
!#  This section contains the initialization subroutines for the info object
!#  as well as the subroutines related to the object itself.
!############################################################################

!### Modules to be Included###
USE Mod_InfoDefs
	
IMPLICIT NONE

CONTAINS

!### Initializes the Info Object ###
SUBROUTINE Info_Init(info)
	Type(Infoc), INTENT(INOUT) :: info

	!### Read in Namelist ###
	CALL Info_AuraData(info)

	!### Initialize Physical Units ###
	CALL Phys_Init(info%physo)

	!### Sanity Checks ###
	if (info%cro%nbins .LT. 2) then
		write(*,*) 'INFO ERROR: Number of CR bins cannot be less than 2'
		STOP
	endif

	if (info%envo%nphbins .NE. 0 .AND. info%envo%nphbins .LT. 2) then
		write(*,*) 'INFO ERROR: Number of ambient photon field bins cannot be less than 2'
		STOP
	endif	

	!### Allocate CR Memory ###
	CALL Info_InitCR(info%cro)

	!### Allocate Ambient Photon Field if needed ###
	if (info%envo%nphbins .GT. 2) then
		CALL Info_InitEnv(info%envo)
	endif

END SUBROUTINE Info_Init

!### Clean up Info Object ###
SUBROUTINE Info_Destroy(info)
	Type(Infoc) :: info

	!### De-allocate CR Memory
	CALL Info_DestroyCR(info%cro)

	!### De-allocate Environment Arrays if needed ###
	if (info%envo%nphbins .GT. 2) then
		CALL Info_DestroyEnv(info%envo)
	endif

END SUBROUTINE Info_Destroy

!### Read in Info Namelist ###
SUBROUTINE Info_AuraData(info)
	Type(Infoc), INTENT(INOUT) :: info

	INTEGER :: ncrbins, nphbins
	LOGICAL namelist_exists, verbose
	CHARACTER*20 oroot, odump_path

	!### Definte the namelist ###
	NAMELIST /Information/ ncrbins, nphbins, verbose, oroot, odump_path

	!### Make sure the namelist file exists ###
	INQUIRE(FILE='AURA.data', EXIST=namelist_exists)
	IF (.NOT. namelist_exists) THEN
		PRINT '("INFO ERROR: setup namelist file AURA.data is missing")'
		STOP
	END IF

	!### read in the namelist ###
	OPEN(UNIT=1, FILE='AURA.data', STATUS="OLD", FORM="FORMATTED")
	READ(1, NML=Information)
	CLOSE(1)

	info%cro%nbins=ncrbins
	info%envo%nphbins=nphbins
	info%verbose=verbose
	info%oroot=oroot
	info%odump_path=odump_path

END SUBROUTINE Info_AuraData

!### Allocate Space for CRs ###
SUBROUTINE Info_InitCR(cro)
	TYPE(CR), INTENT(INOUT) :: cro

	ALLOCATE(cro%n(cro%nbins),cro%e(cro%nbins),cro%g(cro%nbins),cro%ig(cro%nbins),cro%g2(cro%nbins),cro%ig2(cro%nbins),cro%te(cro%nbins),cro%ite(cro%nbins))

END SUBROUTINE Info_InitCR

!### De-Allocate Space for CRs ###
SUBROUTINE Info_DestroyCR(cro)
	TYPE(CR), INTENT(INOUT) :: cro

	DEALLOCATE(cro%n,cro%e,cro%g,cro%ig,cro%g2,cro%ig2,cro%te,cro%ite)

END SUBROUTINE Info_DestroyCR

!### Allocate Space for Ambient Photon Field ###
SUBROUTINE Info_InitEnv(envo)
	TYPE(Env), INTENT(INOUT) :: envo

	ALLOCATE(envo%nph(envo%nphbins),envo%nuph(envo%nphbins))

END SUBROUTINE Info_InitEnv

!### De-Allocate Space for Ambient Photon Field ###
SUBROUTINE Info_DestroyEnv(envo)
	TYPE(Env), INTENT(INOUT) :: envo

	DEALLOCATE(envo%nph,envo%nuph)

END SUBROUTINE Info_DestroyEnv

!### Calculates some useful quantities from the CR spectrum ###
SUBROUTINE Info_CRCalc(physo,cro)
	TYPE(CR), INTENT(INOUT) :: cro
	TYPE(Phys), INTENT(IN) :: physo

	REAL*8 :: imc2

	imc2=1d0/(cro%m*physo%c**2)

	!### Calculating the Lorentz Factors for the distribution ###
	cro%g=cro%e*imc2+1d0
	cro%ig=1d0/cro%g

	!### Calculating the Total Energy for the distribution ###
	cro%te=cro%g*cro%m*physo%c**2
	cro%ite=1d0/cro%te

	!### Other useful quantities ###
	cro%g2=cro%g**2
	cro%ig2=1d0/cro%g2

END SUBROUTINE Info_CRCalc

END MODULE Mod_Info
