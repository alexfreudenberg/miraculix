
/* 
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Collection of system specific auxiliary functions

 Copyright (C) 2001 -- 2017 Martin Schlather, 

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
//#include <unistd.h>
//#include "RandomFieldsUtils.h"
//#include "win_linux_aux.h"
//#include "General_utils.h"
//#include "intrinsics.h"
//#include "Solve.h"
#include "utils.h"
#include "MX.h"
#include "xport_import.h"
#include <Basic_utils.h>
#include <General_utils.h>
#include "kleinkram.h"
//#include "own.h"


int *ToIntDummy = NULL;
int  ToIntN = 0;
int *ToIntI(SEXP X, bool *create, bool round) {
  
  if (TYPEOF(X) == INTSXP) {
    *create = false;
    return INTEGER(X);
  }
  if (TYPEOF(X) == LGLSXP) {
    *create = false;
    return LOGICAL(X);
  }
  int len = length(X);
  if (len > 100 || PL > 1)
    UTILSINFO("Better use 'integer' as storage mode (for one of the arguments).");
  int *y;
  if (*create || ToIntN < len) {
    y = (int *) MALLOC(sizeof(int) * len);    
    if (y == NULL) ERR1("not enough memory for an %d vector of integers", len);
    if (!*create) {
      FREE(ToIntDummy);
      ToIntDummy = y;
      ToIntN = len;
    }
  } else y = ToIntDummy;
  double *x = (double *) REAL(X);
  if (round) for (int i=0; i<len; i++) y[i] = (int) ROUND(x[i]);
  else for (int i=0; i<len; i++) y[i] = (int) x[i];
  return y;
}


void freeGlobals() {
  FREE(ToIntDummy);
}

