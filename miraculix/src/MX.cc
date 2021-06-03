
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2018 -- 2019  Martin Schlather

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


#include <R.h>
#include <Rinternals.h>
#include <inttypes.h> 
//#include <Rdefines.h>
#include "AutoMiraculix.h"
#include "MX.h"
#include <General_utils.h>
#include <zzz_RandomFieldsUtils.h>
#include "options.h"
#include "xport_import.h"



bool debugging = false;


Uint Inti(SEXP X, Uint i) {
  switch(TYPEOF(X)) {
  case INTSXP : return (Uint) INTEGER(X)[i];
  case LGLSXP : return (Uint) LOGICAL(X)[i];
  case REALSXP : return (Uint) REAL(X)[i];
  default : ERR("not of numerical type");
  }
}

Rint * GetAddress(Uint *info, Uint where) {
  addr_Uint addalign;
  for (Uint addr=0; addr<UnitsPerAddress; addr++)      
    addalign.u[addr] = info[where + addr];
  return addalign.a;
}

Uint *GetInfoUnchecked(SEXP Code) { // has negative integers included!!
  SEXP Infos = getAttrib(Code, Information);
  if (TYPEOF(Infos) != INTSXP) ERR("Obsolete storage mode");
  return (Uint *) INTEGER(Infos);
}
