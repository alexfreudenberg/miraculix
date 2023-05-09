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


#include "Basic_miraculix.h"
#include "xport_import.h"




#if defined compatibility_to_R_h
#include "options.h"
#include "miraculix.h"
#include "kleinkram.h"
#include "extern.h"
#include "utils_miraculix.h"
#include "errors_messages.h"


#define importfrom "RandomFieldsUtils"

bool ToFalse[1] = { false };

#define PLoffset -10


KEY_type *PIDKEY_M[PIDMODULUS];

void option_NULL(KEY_type *KT, bool keep_messages) {
  //  printf("%d %d\n", geneticsN, messagesN);
  assert(geneticsN==8 && tuningN == 14 && messagesN == 1);
  messages_options m;
  // printf("NULL!\n");
  if (keep_messages)
    MEMCOPY(&m, &(KT->global.messages), sizeof(messages_options));

  MEMCOPY(&(KT->global), &OPTIONS_MIRACULIX, sizeof(option_type));
  // pointer auf NULL setzten
  if (keep_messages)
    MEMCOPY(&(KT->global.messages), &m, sizeof(messages_options));

  int n = KT->global.genetics.ncentered;
  if (n > 0) {
    Long bytes = n * sizeof(*(KT->global.genetics.pcentered));
    KT->global.genetics.pcentered = (double*) MALLOC(bytes);
    MEMCOPY(KT->global.genetics.pcentered, &OPTIONS_MIRACULIX, bytes);
  }
  
  load_utilsoptions(&(KT->global_utils), false);
}


void option_DELETE(KEY_type *KT) {
   // pointer loeschen
  genetics_options *gp = &(KT->global.genetics);
  FREE(gp->pcentered);
  Ext_del_utilsoption(&(KT->global_utils));
}


void KEY_type_M_NULL(KEY_type *KT) {
  // ACHTUNG!! setzt nur die uninteressanten zurueck. Hier also gar ncihts.
  KT->next = NULL; // braucht es eigentlich nicht
  KT->doshow = true;
  KT->ToIntDummy = NULL;
  KT->ToIntN = 0;
  option_NULL(KT, false);
}

void FREEglobal() {
  KEY_type *KT = KEYT_M(); 
  FREE(KT->ToIntDummy);
  KT->ToIntN = 0;
}

void KEY_type_M_DELETE(KEY_type **S) {
  KEY_type *KT = *S;
  option_DELETE(KT);
  FREEglobal();
  UNCONDFREE(*S);
}

KEY_type *KEYT_M() {
  int mypid;
  Ext_pid(&mypid);
   KEY_type *p = PIDKEY_M[mypid % PIDMODULUS];
  //   printf("%d %d %ld\n", mypid,  PIDMODULUS, p);
  if (p == NULL) {
     KEY_type *neu = (KEY_type *) XCALLOC(1, sizeof(KEY_type));
     assert(neu != NULL);

    PIDKEY_M[mypid % PIDMODULUS] = neu;
    neu->visitingpid = mypid;    
    if (PIDKEY_M[mypid % PIDMODULUS] != neu) { // another process had same idea
      FREE(neu);
      return KEYT_M(); // ... and try again
    }
    neu->pid = mypid;
    neu->visitingpid = 0;
    neu->ok = true;
    if (PIDKEY_M[mypid % PIDMODULUS] != neu) BUG;
    KEY_type_M_NULL(neu);       
    return neu;
  }

   
  while (p->pid != mypid && p->next != NULL) {
      //    printf("pp = %d\n", p->pid);
    p = p->next;
  }
  //  printf("pp m = %d %d\n", p->pid, mypid);
  if (p->pid != mypid) {
    if (!p->ok || p->visitingpid != 0) {
      PRINTF("pid collision %d %d\n",  p->ok, p->visitingpid);
      BUG;
      return KEYT_M();
    }
    p->visitingpid = mypid;
    p->ok = false;
    if (p->visitingpid != mypid || p->ok) {
      return KEYT_M();
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
    KEY_type_M_NULL(neu); 
   return KEYT_M();
  }
   return p;
}

SEXP setlocalRFutils(SEXP seed, SEXP printlevel) {
  KEY_type *KT = KEYT_M();
  assert(KT != NULL);
  utilsoption_type *global_utils = &(KT->global_utils);
  basic_options *basic = &(global_utils->basic);
  assert(global_utils != NULL);
  if (LENGTH(seed) > 0) basic->seed = Integer(seed, (char *) "seed", 0);
  if (LENGTH(printlevel) > 0) {
    basic->Rprintlevel = Integer(printlevel, (char *) "printlevel", 0);
    basic->Cprintlevel = basic->Rprintlevel + PLoffset;
  }
  return R_NilValue;
}


#ifdef __cplusplus
extern "C" {
#endif 
SEXP copyoptions() {
  KEY_type *KT = KEYT_M();
  option_NULL(KT, true);
  return R_NilValue;
}  
#ifdef __cplusplus
}
#endif



void finalizeoptions(int local) {
  if (!local) {
    option_type *options;
    utilsoption_type *utils;
   WhichOptionList(local, &options, &utils);
  }
}

void setoptionsRFU(int i, int j, SEXP el, char name[LEN_OPTIONNAME], 
		       bool isList, bool local) {
  option_type *options;
  utilsoption_type *utils;
  WhichOptionList(local, &options, &utils);
 if (!local && Ext_parallel()) 
    ERR1("Option '%.25s' can be set only through 'RFoptions' at global level",
	 allMoptions[i][j]);
  Ext_setoptionsRFU(i, j, el, name, isList, utils);
}


void getoptionsRFU(SEXP sublist, int i, bool local) {  
  //printf("get local RFU %d\n", i);
  option_type *options;
   utilsoption_type *utils;
   WhichOptionList(local, &options, &utils);
   Ext_getoptionsRFU(sublist, i, utils);
}


// RFU INCLUDE extern "C" void loadoptionsRFU();

void startMiraculix(int x);

  
void loadoptions(int *n) {
#if ! defined DO_PARALLEL
  startMiraculix(1);
#else
  startMiraculix(*n != NA_INTEGER && *n > 0 ? *n : 1);
#endif
  install_default();
    
  MEMSET(PIDKEY_M, 0, PIDMODULUS * sizeof(KEY_type *));
   
  Ext_attachRFUoptions((char *) "miraculix", prefixMlist, prefixMN,
		       allMoptions, allMoptionsN,
		       setMoptions, getMoptions, finalizeoptions, NULL,
		       setoptionsRFU, getoptionsRFU,
		       -10, false,
		       GPU_NEEDS, // from configure.ac
		       check_intrinsics(),
		       MIRACULIX_VERSION,
		       RFU_VERSION,
		       MEMisALIGNED
		       );
  // Ext_attachSetNGet((char*) "miraculix", (char *) "RandomFieldsUtils",
  //		    setoptionsRFU, getoptionsRFU);
  finalizeoptions(false);
}


void PIDKEY_M_DELETE() {
  for (int kn=0; kn<PIDMODULUS; kn++) {
    KEY_type *KT = PIDKEY_M[kn];
    while (KT != NULL) {
      KEY_type *q = KT;
      KT = KT->next;
      KEY_type_M_DELETE(&q);
    }
    PIDKEY_M[kn] = NULL;
  }
}


void detachoptions() {
  PIDKEY_M_DELETE();
  Ext_detachRFUoptions(prefixMlist, prefixMN);
}

SEXP MiraculixOptions(SEXP options) {
  return Ext_RFUoptions(options, (char*) "miraculix");
}


#endif
