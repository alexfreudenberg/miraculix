/*
 Authors 
 Martin Schlather, martin.schlather@uni-mannheim.de

 Copyright (C) 2022-2023 Martin Schlather 

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
*/


#include "Basic_RFUlocal.h" // must be before anything else

#if defined compatibility_to_R_h
#include "RandomFieldsUtils.h"
#include "zzz_RFU.h"

void sortingFromTo(double *d, int len, int from, int to, usr_bool NAlast);
void sortingIntFromTo(int *d, int len, int from, int to,  usr_bool NAlast);


SEXP sortX(SEXP Data, SEXP From, SEXP To, SEXP NAlast) {
  if (length(Data) > MAXINT) BUG;
  int 
    err = NOERROR,
    len = length(Data),
    from = MAX(1, INTEGER(From)[0]),
    to = MIN(INTEGER(To)[0], len);
  if (from > to) return R_NilValue; 

  usr_bool nalast;
  if (LOGICAL(NAlast)[0] == NA_LOGICAL) nalast = Nan;
  else nalast = LOGICAL(NAlast)[0] ? True : False;
  SEXP Ans;

  if (TYPEOF(Data) == REALSXP) {
    //     printf("%d %d %d %d\n",  from, to, INTEGER(To)[0], len);
    PROTECT(Ans=allocVector(REALSXP, to - from + 1));
    double *data;
    int bytes = len * sizeof(*data);
    if ((data = (double*) MALLOC(bytes)) == NULL) { 
      err = ERRORMEMORYALLOCATION; goto ErrorHandling; 
    }
    MEMCOPY(data, REAL(Data), bytes);
    sortingFromTo(data, len, from, to, nalast);
    from--;
    double *ans;
    ans = REAL(Ans);
    for (int i=from; i<to; i++) ans[i - from] = data[i];
    FREE(data);
  } else if  (TYPEOF(Data) == INTSXP) {
    PROTECT(Ans=allocVector(INTSXP, to - from + 1));
    int *data,
      bytes = len * sizeof(*data);
    if ((data = (int*) MALLOC(bytes)) == NULL) { 
      err = ERRORMEMORYALLOCATION; goto ErrorHandling; 
    }
    MEMCOPY(data, INTEGER(Data), bytes);
    sortingIntFromTo(data, len, from, to, nalast);
    from--;
    int *ans;
    ans = INTEGER(Ans);
    for (int i=from ; i<to; i++) ans[i - from] = data[i];
    FREE(data);
  } else ERR0("Data must be real valued or integer valued.");

  
 ErrorHandling :
  UNPROTECT(1);

  switch(err) {
  case ERRORMEMORYALLOCATION : ERR0("not enough memory");
  default:;
  }

  return Ans;
}
 

void orderingFromTo(double *d, int len, int dim, int *pos, int from, int to,
		    usr_bool NAlast);
void orderingIntFromTo(int *d, int len, int dim, int *pos, int from, int to, 
		       usr_bool NAlast);
SEXP orderX(SEXP Data, SEXP From, SEXP To, SEXP NAlast) {
  if (length(Data) > MAXINT) BUG;
  int 
    err = NOERROR,
    len = length(Data),
    from = MAX(1, INTEGER(From)[0]),
    to = MIN(INTEGER(To)[0], len);
  if (from > to) return R_NilValue; 

  SEXP Ans;
  PROTECT(Ans=allocVector(INTSXP, to - from + 1));

  usr_bool nalast;
  if ( LOGICAL(NAlast)[0] == NA_LOGICAL) nalast = Nan;
  else nalast = LOGICAL(NAlast)[0] ? True : False;
  int *pos = (int*) MALLOC(len * sizeof(int));
  if (pos == NULL) {err = ERRORMEMORYALLOCATION; goto ErrorHandling;}

  if (TYPEOF(Data) == REALSXP) {
    //     printf("%d %d %d %d\n",  from, to, INTEGER(To)[0], len);
    orderingFromTo(REAL(Data), len, 1, pos, from, to, nalast);
  } else if  (TYPEOF(Data) == INTSXP) {
    orderingIntFromTo(INTEGER(Data), len, 1, pos, from, to, nalast);
  } else {
    err = ERRORFAILED;
    goto ErrorHandling;
  }

  from--;
  
  int *ans;
  ans = INTEGER(Ans);
  for (int i=from; i<to; i++) ans[i - from] = pos[i] + 1;

 ErrorHandling :
  FREE(pos);
  UNPROTECT(1);
  
  switch(err) {
  case ERRORFAILED : ERR0("Data must be real valued or integer valued.");
  case ERRORMEMORYALLOCATION : ERR0("not enough memory");
  default:;
  }

  return Ans;
}
 


/* 
   extendable to higher dim :

  if (from > to) return R_NilValue; 

  int
    *pos = (int*) MALLOC(len * sizeof(int));
  usr_bool nalast = LOGICAL(NAlast)[0] == NA_LOGICAL ? Nan :
    LOGICAL(NAlast)[0] ? True : False;
  SEXP Ans;


  if (TYPEOF(Data) == REALSXP) {
    // printf("%d %d %d %d\n",  from, to, INTEGER(To)[0], len);
    PROTECT(Ans=allocVector(REALSXP, to - from + 1));
    double *ans = REAL(Ans),
      *data = REAL(Data);
    ordering(data, len, dim, pos, from, to, nalast);
    from--;
    for (int i=from; i<to; i++) {
      //printf("%d %d %d %10g     ", i, from, pos[i], data[pos[i]]);
      ans[i - from] = data[pos[i]];
    }
  } else if  (TYPEOF(Data) == INTSXP) {
    PROTECT(Ans=allocVector(INTSXP, to - from + 1));
    int *ans = INTEGER(Ans),
      *data = INTEGER(Data);
    orderingInt(data, len, dim, pos, from, to, nalast);
    from--;
    for (int i=from ; i<to; i++) ans[i - from] = data[pos[i]];
  } else ERR0("Data must be real valued or integer valued.");
  UNPROTECT(1);
  f ree(pos);
  return Ans;
}
 
*/

#endif
