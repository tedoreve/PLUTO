#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_Thread_H_
#define _CH_Thread_H_

#include "BaseNamespaceHeader.H"

#ifdef _OPENMP
#include <omp.h>
#endif

extern bool onThread0();

#include "BaseNamespaceFooter.H"

#endif
