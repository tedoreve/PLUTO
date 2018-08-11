/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Stellar wind test problem.

  Sets initial condition for a spherically symmetric radial wind blowing
  from the origin of coordinates, see [Mig14]
  The initial condition consists of a constant-density ambient medium with
  \f[
     \rho = \rho_a \,,\qquad  p = p_a\,,\qquad \vec{v} = v_{csm}\hvec{j}
  \f]
  where \c v_csm is the velocity of the star with respect to the background.
  The wind is injeted using the \c INTERNAL_BOUNDARY where flow quantities are
  kept constant in time and equal to
  \f[
     r^2v_r\rho = \rho_0 V_0^2r_0^2 \,,\qquad
     v_r = V_0\hvec{r} \,,\qquad
     p   = \frac{c_s^2}{\Gamma}\rho^\Gamma
  \f]
  These value are defined through the  UserDefBoundary() function when
  \c side is equal to 0.

  Dimensions are chosen so that the spherical wind shell has radius 1,
  density 1 and velocity 1 (\f$ r_0 = V_0 = \rho_0 = 1\f$).

  The input parameters that control the problem dynamics are

  - <tt>g_inputParam[CS_WIND]</tt>: sets the sound speed in the wind region;
  - <tt>g_inputParam[RHO_AMB]</tt>: sets the ambient density
  - <tt>g_inputParam[CS_AMB]</tt>:  sets the ambient sound speed;
  - <tt>g_inputParam[V_CSM]</tt>:   sets the velocity of the star with
                                    respect to the background;

  Configurations #01-05 and #07-08 work in 2D cylindrical axisymmetric
  coordinates while conf. #06 is 3D Cartesian.
  An AMR setup is available with configuration #04.

  \image html hd_stellar_wind.08.jpg "Density map at the end of computation for configuration #8"

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 25, 2012

  \b References
     - [Mig14]: "High-order conservative reconstruction schemes for finite volume methods in cylindrical and spherical coordinates",
       Mignone, JCP (2014) 270, 784
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{
  double rho, cs, R;
  static int first_call = 1;
  

  
  rho = g_inputParam[RHO_AMB];
  cs  = g_inputParam[CS_AMB];
  us[RHO] = rho;
  us[PRS] = cs*cs*rho/g_gamma;
  


  #if GEOMETRY == CARTESIAN
   us[VX1] = 0.0;
   us[VX2] = -g_inputParam[V_AMB];
   //us[VX3] = -g_inputParam[V_AMB];
   us[BX1] = g_inputParam[MAG1];
   us[BX2] = g_inputParam[MAG2];
   //us[BX3] = g_inputParam[MAG3];
  #elif GEOMETRY == CYLINDRICAL
   us[VX1] = 0.0;
   us[VX2] = -g_inputParam[V_AMB];
   us[BX1] = g_inputParam[MAG1];
   us[BX2] = g_inputParam[MAG2];
  #elif GEOMETRY == SPHERICAL
   us[VX1] =  g_inputParam[V_AMB]*cos(x2);
   us[VX2] = -g_inputParam[V_AMB]*sin(x2);
   us[BX1] = g_inputParam[MAG1];
   us[BX2] = g_inputParam[MAG2];
  #endif
  
  #if ADD_BACKGROUND == YES
  if (first_call){
    int k, input_var[256];
    for (k = 0; k< 256; k++) input_var[k] = -1;
    input_var[0] = RHO;    
    input_var[1] = VX1;
    input_var[2] = VX2; 
    input_var[3] = -1;
    InputDataSet ("./grid.out",input_var);
    InputDataRead("./data.0010.dbl"," ");
    first_call = 0;
  }
  InputDataInterpolate(us, x1, x2, x3);  /* -- interpolate density from
                                              input data file -- */
  #endif
  
   if (x2 <= 0.08 && x2 >= 0.03)
  {

    us[RHO] = rho*1.5;
	us[PRS] = cs*cs*rho*1.5/g_gamma;

  }

  R = sqrt(x1*x1 + x2*x2);
  us[TRC] = (R <= 1.0 ? 1.0:0.0);
  g_smallPressure = 1.e-5;
}

#if (BODY_FORCE & POTENTIAL)
/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{   
	double M_star, G;
	G = 1.0e-12;                          //gravitational constant in this unit setting
	M_star  = g_inputParam[M_STAR];       //the mass of the central black hole
	#if GEOMETRY == CARTESIAN
	return -G*M_star/sqrt(x1*x1+x2*x2+x3*x3);
	#elif GEOMETRY == CYLINDRICAL
	return -G*M_star/sqrt(x1*x1+x2*x2);
	#elif GEOMETRY == SPHERICAL
	return -G*M_star/x1;
	#endif
}
#endif

/* **************************************************************** */
void Analysis (const Data *d, Grid *grid)
/*
 *
 ****************************************************************** */
{

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox * box, int side, Grid *grid)
/*
 * Sets inflow boundary condition at the top boundary (side == X2_END)
 * and the stellar wind region when side == 0.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x1, *x2, *x3;
  double  r, r0, cs;
  double  Vwind , rho, vr;

  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;

  Vwind = g_inputParam[V_WIND];
  rho   = g_inputParam[RHO_WIND];
  r0    = g_inputParam[R_WIND];
  cs    = g_inputParam[CS_WIND];
  
  if (side == 0){

    TOT_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN
       r  = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
       if (r <= r0){
         // vr    = tanh(r/r0/0.1)*Vwind;
         // rho   = Vwind*r0*r0/(vr*r*r);
         d->Vc[RHO][k][j][i] = rho;
         d->Vc[VX1][k][j][i] = Vwind*x1[i]/r;
         d->Vc[VX2][k][j][i] = Vwind*x2[j]/r;
         //d->Vc[VX3][k][j][i] = Vwind*x3[k]/r;
         d->Vc[PRS][k][j][i] = cs*cs/g_gamma*rho;
         d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;
       }
      #elif GEOMETRY == CYLINDRICAL
       r  = sqrt(x1[i]*x1[i] + x2[j]*x2[j]);
       if (r <= r0){
         // vr    = tanh(r/r0/0.1)*Vwind;
         // rho   = Vwind*r0*r0/(vr*r*r);
         d->Vc[RHO][k][j][i] = rho;
         d->Vc[VX1][k][j][i] = Vwind*x1[i]/r;
         d->Vc[VX2][k][j][i] = Vwind*x2[j]/r;
         d->Vc[PRS][k][j][i] = cs*cs/g_gamma*rho;
         d->flag[k][j][i]   |= FLAG_INTERNAL_BOUNDARY;
       }
      #endif
    }
  }

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
   X2_BEG_LOOP(k,j,i){ }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */

    cs  = g_inputParam[CS_AMB];
    rho = g_inputParam[RHO_AMB];
    X2_END_LOOP(k,j,i){
      #if GEOMETRY == CYLINDRICAL
       d->Vc[VX1][k][j][i] = 0.0;
       d->Vc[VX2][k][j][i] = -g_inputParam[V_CSM];
       d->Vc[RHO][k][j][i] =  rho;
       d->Vc[PRS][k][j][i] =  cs*cs*rho/g_gamma;
      #endif
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    X3_BEG_LOOP(k,j,i){}
  }

  if (side == X3_END){  /* -- X3_END boundary -- */

    cs  = g_inputParam[CS_AMB];
    rho = g_inputParam[RHO_AMB];
    X3_END_LOOP(k,j,i){
      #if GEOMETRY == CARTESIAN
       d->Vc[VX1][k][j][i] = 0.0;
       d->Vc[VX2][k][j][i] = 0.0;
       d->Vc[VX3][k][j][i] = -g_inputParam[V_AMB];
       d->Vc[RHO][k][j][i] =  rho;
       d->Vc[PRS][k][j][i] =  cs*cs*rho/g_gamma;
      #endif
    }
  }
}
