PROGRAM aura
!############################################################################
!#  This is the program driver for Aura.
!############################################################################

!### Modules to be Included ###
USE Mod_Info
USE Mod_InfoDefs
USE Mod_Emis
USE Mod_User
USE Mod_IO

IMPLICIT NONE

	!### Declare our objects ###
	TYPE(Infoc) :: info
	TYPE(Emis) :: emiso

	REAL :: tstart, tend

	CALL CPU_TIME(tstart)

	write(*,*) 'Starting Aura'

	!### Initializing Info ###
	CALL Info_Init(info)

	!### Initializing Emis ###
	CALL Emis_Init(emiso,info)

	!### Define Initial Conditions ### 
	CALL User_Init(info)

	!### Calculate Emissions ###
	CALL Emis_Calc(info,emiso)

	!### Output Results ###
	CALL IO_Emis(emiso,info)

	!### Destroy Emis ###
	CALL Emis_Destroy(emiso)

	!### Destroy Info ###
	CALL Info_Destroy(info)

	CALL CPU_TIME(tend)

	write(*,*) 'Calculation took:', tend-tstart, 'seconds'

END PROGRAM aura
