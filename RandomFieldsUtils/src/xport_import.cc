/*
 Authors
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2015 -- 2019 Martin Schlather, 

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

#include <R_ext/Rdynload.h>
#include "intrinsics.h"
#include "options.h"
#include "xport_import.h"
#include "General_utils.h"
#include "win_linux_aux.h"
#include "zzz_RandomFieldsUtils.h"
#include "RandomFieldsUtils.h"
//#include "zzz_RandomFieldsUtils.h" // always last


#define PLoffset -10
int PL = C_PRINTLEVEL,
  CORES = 1; //  return;  TO DO: replace by KEYT->global_utils


KEY_type *PIDKEY[PIDMODULUS];
int parentpid=0;
bool parallel() {
  int mypid;
  pid(&mypid);
  return mypid != parentpid;
}


void KEY_type_NULL(KEY_type *KT) {
  // ACHTUNG!! setzt nur die uninteressanten zurueck. Hier also gar ncihts.
  KT->next = NULL; // braucht es eigentlich nicht
  KT->doshow = true;
  KT->ToIntDummy = NULL;
  KT->ToIntN = 0;
  KT->ToRealDummy = NULL;
  KT->ToRealN = 0;
  KT->nu2old = KT->nuOld = KT->nu1old = KT->nuAlt = -RF_INF;
  // option_type_NULL(KT, false)
}

void KEY_type_DELETE(KEY_type **S) {
  KEY_type *KT = *S;
  //option_type_DELETE(KT);
  FREE(KT->ToIntDummy);
  FREE(KT->ToRealDummy);
  UNCONDFREE(*S);
}



KEY_type *KEYT() {
  int mypid;
  pid(&mypid);
  //  printf("entering KEYT %d %d \n", mypid, parentpid);
  // for (int i=0; i<PIDMODULUS; i++) printf("%ld ", PIDKEY[i]);
  KEY_type *p = PIDKEY[mypid % PIDMODULUS];
  //  printf("%d %d %ld\n", mypid,  PIDMODULUS, p);
  if (p == NULL) {
    KEY_type *neu = (KEY_type *) XCALLOC(1, sizeof(KEY_type));
    assert(neu != NULL);
    PIDKEY[mypid % PIDMODULUS] = neu;
    neu->visitingpid = mypid;    
    if (PIDKEY[mypid % PIDMODULUS] != neu) { // another process had the
      //                                        same idea
      FREE(neu);
      return KEYT(); // ... and try again
    }
    neu->pid = mypid;
    //    printf("neu %d %d\n", mypid);
    neu->visitingpid = 0;
    neu->ok = true;
    if (PIDKEY[mypid % PIDMODULUS] != neu) BUG;
    KEY_type_NULL(neu);    
    WARN_PARALLEL;
   return neu;
  }
  while (p->pid != mypid && p->next != NULL) {
    //    printf("pp = %d\n", p->pid);
    p = p->next;
  }
  //  printf("pp m = %d %d\n", p->pid, mypid);
  if (p->pid != mypid) {
    if (!p->ok || p->visitingpid != 0) {
      if (PL >= PL_ERRORS) {
	PRINTF("pid collision %d %d\n",  p->ok, p->visitingpid);
      }
      //
      BUG;
      return KEYT();
    }
    p->visitingpid = mypid;
    p->ok = false;
    if (p->visitingpid != mypid || p->ok) {
      return KEYT();
    }
    KEY_type *neu = (KEY_type *) XCALLOC(1, sizeof(KEY_type));
    neu->pid = mypid;
    if (!p->ok && p->visitingpid == mypid) {
      p->next = neu;
      p->visitingpid = 0;
      p->ok = true;      
      return neu;
    }
    FREE(neu);
    p->visitingpid = 0;
    p->ok = true;
    KEY_type_NULL(neu); 
   return KEYT();
  }
  return p;
}


void setoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], bool isList, bool local);
void getoptions(SEXP sublist, int i, bool local);
void deloptions(bool local);


#define NEED_GPU true
#define NEED_AVX2 true
#define NEED_AVX true
  //#define NEED_SSE4 true
#define NEED_SSSE3 false
#define NEED_SSE2 true
#define NEED_SSE false
void loadoptions() {
  for (int i=0; i<PIDMODULUS; i++) PIDKEY[i] = NULL; 
  pid(&parentpid);
  //  printf("install needd %d\n", NEEDS);
  attachRFoptions((char *) "RandomFieldsUtils",
		  prefixlist, prefixN,
		  all, allN,
  		  setoptions,
		  NULL, // final
		  getoptions,
		  deloptions,
		  0, true, AVX_NEEDS, GPU_NEEDS, AVX_INFO);
  //finalizeoptions();
  SetLaMode();
}


void PIDKEY_DELETE() {
  for (int kn=0; kn<PIDMODULUS; kn++) {
    KEY_type *KT = PIDKEY[kn];
    while (KT != NULL) {
      KEY_type *q = KT;
      KT = KT->next;
      KEY_type_DELETE(&q);
    }
    PIDKEY[kn] = NULL;
  }
}


void  detachoptions(){
  PIDKEY_DELETE();
  detachRFoptions(prefixlist, prefixN);
}


