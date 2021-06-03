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
//#include "RF.h"
#include "IntrinsicsBase.h"
#include <General_utils.h>
#include "options.h"
#include "xport_import.h"
#include "miraculix.h"
#include "mmagpu.h"
#include "kleinkram.h"


#define importfrom "RandomFieldsUtils"

#ifdef CALL
#undef CALL
#endif
#define CALL(what) what##_type Ext_##what = NULL
UTILSCALLS;

#undef CALL
#define CALL(what) Ext_##what = (what##_type) R_GetCCallable(importfrom, #what)
void includeXport() {
  UTILSCALLS;
} // export C




bool ToFalse[1] = { false };

#define PLoffset -10
int PL = C_PRINTLEVEL,
  CORES = INITCORES; //  return;  TO DO: replace by KEYT->global_utils

SEXP Information = NULL,
  Coding = NULL;
bool HAS_CUDA = false;


utilsoption_type *OPTIONS_UTILS;
KEY_type *PIDKEY[PIDMODULUS];
int parentpid=0;
bool parallel() {
  int mypid;
  //  printf("pid\n");
  Ext_pid(&mypid);
  //  printf("pid = %d %d\n", mypid, parentpid);
  return mypid != parentpid;
}



void option_NULL(KEY_type *KT, bool keep_messages) {
  //  printf("%d %d\n", geneticsN, messagesN);
  assert(geneticsN==6 && messagesN == 1);
  messages_options m;
  if (keep_messages)
    MEMCOPY(&m, &(KT->global.messages), sizeof(messages_options));

  MEMCOPY(&(KT->global), &OPTIONS, sizeof(option_type));
  // pointer auf NULL setzten
  if (keep_messages)
    MEMCOPY(&(KT->global.messages), &m, sizeof(messages_options));

  int n = KT->global.genetics.ncentered;
  if (n > 0) {
    int bytes = n * sizeof(double);
    KT->global.genetics.pcentered = (double*) MALLOC(bytes);
    MEMCOPY(KT->global.genetics.pcentered, &OPTIONS, bytes);
  }
  
  Ext_utilsoption_NULL(&(KT->global_utils));
}


void option_DELETE(KEY_type *KT) {
   // pointer loeschen
  genetics_options *gp = &(KT->global.genetics);
  FREE(gp->pcentered);
  FREE(KT->ToIntDummy);
  Ext_utilsoption_DELETE(&(KT->global_utils));
}


void KEY_type_NULL(KEY_type *KT) {
  // ACHTUNG!! setzt nur die uninteressanten zurueck. Hier also gar ncihts.
  KT->next = NULL; // braucht es eigentlich nicht
  KT->doshow = true;
  KT->ToIntDummy = NULL;
  KT->ToIntN = 0;
  option_NULL(KT, false);
}

void KEY_type_DELETE(KEY_type **S) {
  KEY_type *KT = *S;
  option_DELETE(KT);
  FREE(KT->ToIntDummy);
  UNCONDFREE(*S);
}

KEY_type *KEYT() {
  int mypid;
  Ext_pid(&mypid);
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
    if (OPTIONS_UTILS->basic.warn_parallel && mypid == parentpid) {
      PRINTF("Do not forget to run 'RFoptions(storing=FALSE)' after each call of a parallel command (e.g. from packages 'parallel') that calls a function in 'RandomFields'. (OMP within RandomFields is not affected.) This message can be suppressed by 'RFoptions(warn_parallel=FALSE)'.\n"); // ok
    }
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


#ifdef __cplusplus
extern "C" {
#endif 
SEXP copyoptions() {
  KEY_type *KT = KEYT();
  option_NULL(KT, true);
  return R_NilValue;
}

SEXP setlocalRFutils(SEXP seed, SEXP printlevel) {
  KEY_type *KT = KEYT();
  assert(KT != NULL);
  utilsoption_type *global_utils = &(KT->global_utils);
  assert(global_utils != NULL);
  if (length(seed) > 0)
    global_utils->basic.seed = Integer(seed, (char *) "seed", 0);
  if (length(printlevel) > 0) {
    PL = global_utils->basic.Rprintlevel =
      Integer(printlevel, (char *) "printlevel", 0);
    global_utils->basic.Cprintlevel = global_utils->basic.Rprintlevel +PLoffset;
  }
  return R_NilValue;
}
  
#ifdef __cplusplus
}
#endif


void finalizeoptions() {
  utilsoption_type *global_utils = OPTIONS_UTILS;
  PL = global_utils->basic.Cprintlevel - PLoffset;
  CORES = global_utils->basic.cores;
}

void setoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME], bool isList, bool local);
void getoptions(SEXP sublist, int i, bool local);


extern void MultiplyInit128();
extern void MultiplyInit256();
extern void PackedInit128();
extern void PackedInit256();
extern void ShuffleInit128();
extern void ShuffleInit256();
void Init() {
  MultiplyInit128();
  MultiplyInit256();
  PackedInit128();
  PackedInit256();
  ShuffleInit128();
  ShuffleInit256();
}


#define NEED_GPU true
#define NEED_AVX2 true
#define NEED_AVX false
#define NEED_SSSE3  true
#define NEED_SSE2 true
#define NEED_SSE false
void loadoptions() {  
  assert(Haplo == 30);
  for (int i=0; i<PIDMODULUS; i++) PIDKEY[i] = NULL;
  //printf("loading..\n");
  Information = install("information");
  Coding = install("coding");
  HAS_CUDA = check_cuda();

  if (BytesPerUnit != sizeof(Uint))
    ERR1("The programme relies on the assumption that an unsigned integer has %d Bytes.", BytesPerUnit);
  
  if (sizeof(int) != sizeof(Uint))
    ERR2("The programme relies on the assumption that a signed integer the same lenth than an unsigned integer. Found %d and %d bytes respectively.", sizeof(int), sizeof(Uint));
  
  if (sizeof(uint64_t) != sizeof(double))
    ERR2("The programme relies on the assumption that an 'uint64_t' has the size of a 'double'. Found %d and %d bytes respectively.", sizeof(uint64_t), sizeof(double));
  
  if (sizeof(uint64_t) != 8) ERR1("The programme relies on the assumption that an 'uint64_t' has 8 Bytes. Found %d bytes.", sizeof(uint64_t));
  
  if (sizeof(void*) > MaxUnitsPerAddress * BytesPerUnit)
    ERR2("The programme relies on the assumption that an 'void*' has at most %d Bytes. Found %d bytes.", MaxUnitsPerAddress * BytesPerUnit, sizeof(void*));

  includeXport();
  Ext_pid(&parentpid);
  Ext_getUtilsParam(&OPTIONS_UTILS);
  utilsoption_type *global_utils = OPTIONS_UTILS; // OK
  global_utils->solve.max_chol = 8192;
  global_utils->solve.max_svd = 6555; 
  Ext_attachRFoptions((char *) "miraculix", prefixlist, prefixN, all, allN,
		      setoptions, finalizeoptions, getoptions,
		      NULL, -10, false, AVX_NEEDS, GPU_NEEDS, AVX_INFO);
  finalizeoptions();
  Init(); 
}




option_type *WhichOptionList(bool local) {  
  if (local) {
    KEY_type *KT = KEYT();
    if (KT == NULL) BUG;
    //    printf("wol = %ld\n", &(KT->global));
    return &(KT->global);
  }
  return &OPTIONS;
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


void detachoptions() {
  PIDKEY_DELETE();
  Ext_detachRFoptions(prefixlist, prefixN);
}


