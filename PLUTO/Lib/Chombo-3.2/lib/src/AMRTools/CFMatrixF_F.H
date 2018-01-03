#ifndef _CFMATRIXF_F_H_
#define _CFMATRIXF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETCOARSEFINEINTERPMATRIX 
#define GUARDGETCOARSEFINEINTERPMATRIX 
// Prototype for Fortran procedure GETCOARSEFINEINTERPMATRIX ...
//
void FORTRAN_NAME( GETCOARSEFINEINTERPMATRIX ,getcoarsefineinterpmatrix )(
      CHFp_FRA(cfmat)
      ,CHFp_CONST_INT(nref)
      ,CHFp_CONST_INT(nnbrs)
      ,CHFp_CONST_INT(npwrs)
      ,CHFp_CONST_INT(degree)
      ,CHFp_CONST_INT(nfine)
      ,CHFp_CONST_VI(nbrsv)
      ,CHFp_BOX(fbox)
      ,CHFp_BOX(degbox) );

#define FORT_GETCOARSEFINEINTERPMATRIX FORTRAN_NAME( inlineGETCOARSEFINEINTERPMATRIX, inlineGETCOARSEFINEINTERPMATRIX)
#define FORTNT_GETCOARSEFINEINTERPMATRIX FORTRAN_NAME( GETCOARSEFINEINTERPMATRIX, getcoarsefineinterpmatrix)

inline void FORTRAN_NAME(inlineGETCOARSEFINEINTERPMATRIX, inlineGETCOARSEFINEINTERPMATRIX)(
      CHFp_FRA(cfmat)
      ,CHFp_CONST_INT(nref)
      ,CHFp_CONST_INT(nnbrs)
      ,CHFp_CONST_INT(npwrs)
      ,CHFp_CONST_INT(degree)
      ,CHFp_CONST_INT(nfine)
      ,CHFp_CONST_VI(nbrsv)
      ,CHFp_BOX(fbox)
      ,CHFp_BOX(degbox) )
{
 CH_TIMELEAF("FORT_GETCOARSEFINEINTERPMATRIX");
 FORTRAN_NAME( GETCOARSEFINEINTERPMATRIX ,getcoarsefineinterpmatrix )(
      CHFt_FRA(cfmat)
      ,CHFt_CONST_INT(nref)
      ,CHFt_CONST_INT(nnbrs)
      ,CHFt_CONST_INT(npwrs)
      ,CHFt_CONST_INT(degree)
      ,CHFt_CONST_INT(nfine)
      ,CHFt_CONST_VI(nbrsv)
      ,CHFt_BOX(fbox)
      ,CHFt_BOX(degbox) );
}
#endif  // GUARDGETCOARSEFINEINTERPMATRIX 

#ifndef GUARDAPPLYCOARSEFINEINTERP 
#define GUARDAPPLYCOARSEFINEINTERP 
// Prototype for Fortran procedure APPLYCOARSEFINEINTERP ...
//
void FORTRAN_NAME( APPLYCOARSEFINEINTERP ,applycoarsefineinterp )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_CONST_FRA(stmat)
      ,CHFp_CONST_INTVECT(fbase)
      ,CHFp_CONST_INTVECT(cbase)
      ,CHFp_CONST_VI(nbrsv)
      ,CHFp_BOX(fbox) );

#define FORT_APPLYCOARSEFINEINTERP FORTRAN_NAME( inlineAPPLYCOARSEFINEINTERP, inlineAPPLYCOARSEFINEINTERP)
#define FORTNT_APPLYCOARSEFINEINTERP FORTRAN_NAME( APPLYCOARSEFINEINTERP, applycoarsefineinterp)

inline void FORTRAN_NAME(inlineAPPLYCOARSEFINEINTERP, inlineAPPLYCOARSEFINEINTERP)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_CONST_FRA(stmat)
      ,CHFp_CONST_INTVECT(fbase)
      ,CHFp_CONST_INTVECT(cbase)
      ,CHFp_CONST_VI(nbrsv)
      ,CHFp_BOX(fbox) )
{
 CH_TIMELEAF("FORT_APPLYCOARSEFINEINTERP");
 FORTRAN_NAME( APPLYCOARSEFINEINTERP ,applycoarsefineinterp )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_CONST_FRA(stmat)
      ,CHFt_CONST_INTVECT(fbase)
      ,CHFt_CONST_INTVECT(cbase)
      ,CHFt_CONST_VI(nbrsv)
      ,CHFt_BOX(fbox) );
}
#endif  // GUARDAPPLYCOARSEFINEINTERP 

}

#endif
