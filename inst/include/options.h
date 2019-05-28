


/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2016 Martin Schlather

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


#ifndef miraculix_options_H
#define miraculix_options_H 1

#include <Rinternals.h>
#include "local1.h"

#define relationshipN 6
typedef struct relationship_param {
  bool normalized, returnsigma, per_snp;
  usr_bool centered;
  double digits;
  int method;
  // internal:
  double *pcentered;
  int ncentered;
} relationship_param;

#define relationship_START {			\
    true, false, false,				\
      True,					\
      3.0,					\
      Shuffle,					\
      NULL, 0 /* internal */			\
}


typedef struct globalparam {
  relationship_param relationship;
} globalparam;
extern globalparam GLOBAL;

#define prefixN 1
extern const char * prefixlist[prefixN], **all[prefixN];
extern int allN[prefixN];
void setparameter(int i, int j, SEXP el, char name[200], bool isList, int local);
void getRFoptions(SEXP *sublist);
void finalparameter(int local);


#endif
