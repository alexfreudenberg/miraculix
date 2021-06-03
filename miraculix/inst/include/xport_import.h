/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2015 -- 2019 Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/


#ifndef Miraculixxport_H
#define Miraculixxport_H 1

#include "local3.h"


#define UTILSCALLS				\
  CALL(solve_NULL);				\
  CALL(solve_DELETE);				\
  CALL(getUtilsParam);				\
  CALL(attachRFoptions);			\
  CALL(detachRFoptions);			\
  CALL(ordering);				\
  CALL(orderingInt);				\
  CALL(solvePosDefSp);				\
  CALL(sleepMicro);				\
  CALL(utilsoption_NULL);			\
  CALL(utilsoption_DELETE);			\
  CALL(scalarX);				\
  CALL(pid)
  
#ifdef CALL
#undef CALL
#endif
#define CALL(what) extern what##_type Ext_##what
UTILSCALLS;

void includeXport();
extern utilsoption_type* OPTIONS_UTILS;
extern int CORES;
extern int PL;
void PIDKEY_DELETE();
typedef
struct option_type option_type;
option_type *WhichOptionList(bool local);
bool parallel();


// #undef ToReal
// #define ToReal(X) (TYPEOF(X) == REALSXP ? REAL(X) : Ext_ToRealI(X, ToFalse))
//#define ToRealX(X, ) (TYPEOF(X) == REALSXP ? REAL(X) : Ext_ToRealI(X, ToFalse))

// extern bool ToFalse[1];
// #define ToInt(X) Ext_ToIntI(X, ToFalse, false) -- nicht verwenden, da
// parallele Aufrufe !!
//#define ToIntX(X,C) Ext_ToIntI(X, C, false)



#define FREEint(X) if (create##X) FREE(X##int)

extern bool HAS_CUDA;
extern SEXP Information, Coding; // READONLY !!!!


#define check_cuda check_7_5 // the lowest number for which code exists


#endif
