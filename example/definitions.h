#define  PHYSICS                 MHD
#define  DIMENSIONS              2
#define  COMPONENTS              2
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     12

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD        NO
#define  RESISTIVITY             NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  E_EJ                    0
#define  M_EJ                    1
#define  R_EJ                    2
#define  N_H                     3
#define  U_AM                    4
#define  S_PI                    5
#define  N_PI                    6
#define  ETA                     7
#define  BMAG                    8
#define  THETA                   9
#define  PHI                     10
#define  GAMMA                   11

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            CONST_mp
#define  UNIT_LENGTH             CONST_pc
#define  UNIT_VELOCITY           1.0e9
#define  ADD_TURBULENCE          YES

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          YES
#define  PRINT_TO_FILE             NO
#define  INTERNAL_BOUNDARY         NO
#define  SHOCK_FLATTENING          NO
#define  CHAR_LIMITING             YES
#define  LIMITER                   VANLEER_LIM
#define  CT_EMF_AVERAGE            ARITHMETIC
#define  CT_EN_CORRECTION          YES
#define  ASSIGN_VECTOR_POTENTIAL   YES
#define  UPDATE_VECTOR_POTENTIAL   NO
