MODULE Mod_Emis
!############################################################################
!#  This section contains the initialization subroutines for the emission object
!#  as well as the subroutines related to the object itself.
!############################################################################

!### Modules to be Included###
USE Mod_Sync
USE Mod_IC
USE Mod_Brem
USE Mod_PP
USE Mod_Info
	
IMPLICIT NONE

!### Define Emission Class ###
TYPE :: Emis

	!### Define the objects contained in this class ###
	TYPE(Sync) :: synco
	TYPE(IC) :: ico
	TYPE(Brem) :: bremo
	TYPE(PP) :: ppo

	!### Array to hold the frequency bins for the emissions ###
	REAL*8, ALLOCATABLE, DIMENSION(:) :: nuph
	
	!### Lower and Upper Bounds of Frequency Range, Secondary Energy Range ###
	REAL*8 :: nulo, nuhi, elo, ehi

	!### Number of Frequency Bins, number of Secondary CR bins ###
	INTEGER :: nphbins, ncrbins

END TYPE

CONTAINS

!### Initializes the Emission Object ###
SUBROUTINE Emis_Init(emiso,info)
	TYPE(Emis), INTENT(INOUT) :: emiso
	TYPE(Infoc), INTENT(IN) :: info

	!### Read in Namelist ###
	CALL Emis_AuraData(emiso)

	!### Sanity Checks ###
	if (emiso%nphbins .LT. 1) then
		write(*,*) 'EMIS ERROR: Number of emission bins cannot be less than 1'
		STOP
	endif

	if (emiso%nulo .GE. emiso%nuhi .AND. emiso%nphbins .GE. 2) then
		write(*,*) 'EMIS ERROR: The lower frequency bound must be less than the upper bound'
		STOP
	endif

	if (emiso%nphbins .EQ. 1 .AND. emiso%nulo .NE. emiso%nuhi) then
		write(*,*) 'EMIS WARNING: For a single frequency bin, the lower bound is used'
	endif

	!### Initialize Frequency Range ###
	CALL Emis_Freq(emiso)

	!### Initialize Synchrotron Object ###
	if (emiso%synco%synchrotron) then
		CALL Sync_Init(emiso%synco,emiso%nuph,emiso%nphbins)
	endif

	!### Initialize IC Object ###
	if (emiso%ico%ic) then
		CALL IC_Init(emiso%ico,emiso%nuph,emiso%nphbins,info%envo%nuph,info%envo%nph,info%envo%nphbins)
	endif

	!### Initialize Bremsstrahlung Object ###
	if (emiso%bremo%bremsstrahlung) then
		CALL Brem_Init(emiso%bremo,emiso%nuph,emiso%nphbins)
	endif

	!### Initialize Proton-Proton Interaction Object ###
	if (emiso%ppo%pp) then
		CALL PP_Init(emiso%ppo,info%physo,emiso%nuph,emiso%nphbins,emiso%elo,emiso%ehi,emiso%ncrbins)
	endif

END SUBROUTINE Emis_Init

!### Clean up Emission Object ###
SUBROUTINE Emis_Destroy(emiso)
	TYPE (Emis), INTENT(INOUT) :: emiso

	!### Clean-up Synchrotron Object ###
	if (emiso%synco%synchrotron) then
		CALL Sync_Destroy(emiso%synco)
	endif

	!### Clean-up IC Object ###
	if (emiso%ico%ic) then
		CALL IC_Destroy(emiso%ico)
	endif

	!### Clean-up Bremsstrahlung Object ###
	if (emiso%bremo%bremsstrahlung) then
		CALL Brem_Destroy(emiso%bremo)
	endif

	!### Clean-up Proton-Proton Interaction Object ###
	if (emiso%ppo%pp) then
		CALL PP_Destroy(emiso%ppo)
	endif

	!### De-Allocate Arrays ###
	DEALLOCATE(emiso%nuph)

END SUBROUTINE Emis_Destroy

!### Read in Emis Namelist ###
SUBROUTINE Emis_AuraData(emiso)
	TYPE(Emis), INTENT(INOUT) :: emiso

	REAL*8 :: nulo, nuhi, elo, ehi
	INTEGER :: nphbins, ncrbins
	LOGICAL namelist_exists, synchrotron, ic, bremsstrahlung, pp

	!### Definte the namelist ###
	NAMELIST /Emission/ nphbins, nulo, nuhi, ncrbins, elo, ehi, synchrotron, ic, bremsstrahlung, pp

	!### Make sure the namelist file exists ###
	INQUIRE(FILE='AURA.data', EXIST=namelist_exists)
	IF (.NOT. namelist_exists) THEN
		PRINT '("EMIS ERROR: setup namelist file AURA.data is missing")'
		STOP
	END IF

	!### read in the namelist ###
	OPEN(UNIT=1, FILE='AURA.data', STATUS="OLD", FORM="FORMATTED")
	READ(1, NML=Emission)
	CLOSE(1)

	emiso%nphbins=nphbins
	emiso%nulo=nulo
	emiso%nuhi=nuhi
	emiso%ncrbins=ncrbins
	emiso%elo=elo
	emiso%ehi=ehi
	emiso%synco%synchrotron=synchrotron
	emiso%ico%ic=ic
	emiso%bremo%bremsstrahlung=bremsstrahlung
	emiso%ppo%pp=pp

END SUBROUTINE Emis_AuraData

!### Sets the frequency range for emissions ###
SUBROUTINE Emis_Freq(emiso)
	TYPE(Emis), INTENT(INOUT) :: emiso

	REAL*8 :: dlognu
	INTEGER :: i

	!### Allocate space for frequency ###
	ALLOCATE(emiso%nuph(emiso%nphbins))

	!### Check to see if there is a single bin ###
	if (emiso%nphbins .EQ. 1) then
		emiso%nuph(1)=emiso%nulo
		RETURN
	endif

	!### Set size of logarithmic spacing ###
	dlognu=(dlog(emiso%nuhi)-dlog(emiso%nulo))/REAL(emiso%nphbins-1,8)

	!### Calculate frequencies in log space ###
	do i=1,emiso%nphbins
		emiso%nuph(i)=dlog(emiso%nulo)+REAL(i-1,8)*dlognu
	enddo

	!### Convert back to normal space ###
	emiso%nuph=dexp(emiso%nuph)

END SUBROUTINE Emis_Freq

!### Runs the emissions calculation ###
SUBROUTINE Emis_Calc(info,emiso)
	TYPE(Infoc), INTENT(INOUT) :: info
	TYPE(Emis), INTENT(INOUT) :: emiso

	if (info%verbose) write(*,*) 'Starting Emissions'

	!### Calculate some initial quantities that will be used repeated for the CR spectrum ###
	CALL Info_CRCalc(info%physo,info%cro)

	!### Synchrotron ###
	if (emiso%synco%synchrotron) then
		if (info%verbose) write(*,*) 'Starting Synchrotron'

		CALL Sync_Calc(emiso%synco,info%cro,info%physo,info%envo%ne,info%envo%B,info%envo%theta)
	endif

	!### IC ###
	if (emiso%ico%ic) then
		if (info%verbose) write(*,*) 'Starting IC' 

		CALL IC_Calc(emiso%ico,info%cro,info%physo)
	endif

	!### Bremsstrahlung ###
	if (emiso%bremo%bremsstrahlung) then
		if (info%verbose) write(*,*) 'Starting Bremsstrahlung' 

		CALL Brem_Calc(emiso%bremo,info%cro,info%physo,info%envo%ne,info%envo%np,info%envo%nalpha,info%envo%nH,info%envo%nHe)
	endif

	!### Proton-Proton Interactions ###
	if (emiso%ppo%pp) then
		if (info%verbose) write(*,*) 'Starting Proton-Proton Interactions' 

		CALL PP_Calc(emiso%ppo,info%cro,info%physo,info%envo%np,info%envo%nH)
	endif

END SUBROUTINE Emis_Calc

END MODULE Mod_Emis
