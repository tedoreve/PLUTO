MODULE Mod_InfoDefs

!############################################################################
!#  This section contains the class definitions for the info object
!############################################################################

!### Modules to be Included###
USE Mod_Phys
	
IMPLICIT NONE

!### Define Cosmic Ray (CR) Class ##
TYPE :: CR

	!### Array to hold CR spectrum, kinetic energy, and lorentz factor bins ###
	REAL*8, ALLOCATABLE, DIMENSION(:) :: n, e, g, ig, g2, ig2, te, ite

	!### Variables to hold what type of CR it is, charge and mass ###
	REAL*8 :: q, m

	!### Hold onto number of bins ###
	INTEGER :: nbins

END TYPE

!### Define Environment Class ###
TYPE :: Env

	!### Array to hold ambient photon spectrum and frequency bins ###
	REAL*8, ALLOCATABLE, DIMENSION(:) :: nph, nuph

	!### Variables to hold on to ambient conditions ###
	REAL*8 :: ne, np, nalpha, nH, nHe, B, theta

	!### Hold on to number of bins ###
	INTEGER :: nphbins

END TYPE

!### Define Info Class ###
TYPE :: Infoc

	!### Define the objects contained in this class ###
	TYPE(Phys) :: physo
	TYPE(CR) :: cro
	TYPE(Env) :: envo

	!### Verbosity Flag ###
	LOGICAL :: verbose=.FALSE.

	!### IO Info ###
	CHARACTER*20 :: oroot, odump_path

END TYPE

END MODULE Mod_InfoDefs
