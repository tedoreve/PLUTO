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
#define  USER_DEF_PARAMETERS     10

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  DIVB_CONTROL            EIGHT_WAVES
#define  BACKGROUND_FIELD        NO
#define  RESISTIVITY             NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  RHO_AMB                 0
#define  CS_AMB                  1
#define  V_AMB                   2
#define  MAG1                    3
#define  MAG2                    4
#define  MAG3                    5
#define  RHO_WIND                6
#define  CS_WIND                 7
#define  V_WIND                  8
#define  R_WIND                  9

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            CONST_mp
#define  UNIT_LENGTH             CONST_pc
#define  UNIT_VELOCITY           1.0e9

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING         NO
#define  WARNING_MESSAGES          NO
#define  PRINT_TO_FILE             NO
#define  INTERNAL_BOUNDARY         YES
#define  SHOCK_FLATTENING          NO
#define  CHAR_LIMITING             NO
#define  LIMITER                   VANALBADA_LIM
#define  ASSIGN_VECTOR_POTENTIAL   NO
#define  UPDATE_VECTOR_POTENTIAL   NO
