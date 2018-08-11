#define  PHYSICS                 MHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     23

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD        NO
#define  RESISTIVITY             NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               EXPLICIT
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  E_EJ                    0
#define  M_EJ                    1
#define  R_EJ                    2
#define  M_STAR                  3
#define  N_H                     4
#define  U_AM                    5
#define  S_PI                    6
#define  N_PI                    7
#define  W_C                     8
#define  BMAG                    9
#define  THETA                   10
#define  PHI                     11
#define  GAMMA                   12
#define  Temp                    13
#define  NU_VISC                 14
#define  RHO_AMB                 15
#define  CS_AMB                  16
#define  V_AMB                   17
#define  DIST                    18
#define  RHO_WIND                19
#define  CS_WIND                 20
#define  V_WIND                  21
#define  R_WIND                  22

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            CONST_mp
#define  UNIT_LENGTH             CONST_pc
#define  UNIT_VELOCITY           1.0e9
#define  ADD_BACKGROUND          NO

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          NO
#define  PRINT_TO_FILE             NO
#define  INTERNAL_BOUNDARY         YES
#define  SHOCK_FLATTENING          NO
#define  CHAR_LIMITING             YES
#define  LIMITER                   VANLEER_LIM
#define  CT_EMF_AVERAGE            ARITHMETIC
#define  CT_EN_CORRECTION          YES
#define  ASSIGN_VECTOR_POTENTIAL   YES
#define  UPDATE_VECTOR_POTENTIAL   NO
