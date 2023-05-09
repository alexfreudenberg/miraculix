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


#include "Basic_RFUlocal.h"

#if defined compatibility_to_R_h

#include <R_ext/Rdynload.h>
#include "options_RFU.h"
#include "zzz_RFU.h"
#include "win_linux_aux.h"
#include "extern_RFU.h"
#include "RandomFieldsUtils.h"
#include "xport_import_RFU.h"
#include "compatibility.SEXP.h"



KEY_type *PIDKEY[PIDMODULUS];

void KEY_type_NULL(KEY_type *KT) {
  // ACHTUNG!! setzt nur die uninteressanten zurueck. Hier also gar ncihts.
  KT->next = NULL; // braucht es eigentlich nicht
  KT->doshow = true;
  KT->ToIntDummy = NULL;
  KT->ToIntN = 0;
  KT->ToRealDummy = NULL;
  KT->ToRealN = 0;
  whittle_work_type *work = &(KT->whittle);
  work->nu2old = work->nuOld = work->nu1old = work->nuAlt = -RF_INF;
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
    
    //if (basic.warn_parallel && mypid == parentpid)  PRINTF("Do not forget to run 'RFoptions(storing=FALSE)' after each call of a parallel command (e.g. from packages 'parallel') that calls a function in 'RandomFields'. (OMP within RandomFields is not affected.) This message can be suppressed by 'RFoptions(warn_parallel=FALSE)'.") /*// OK */
    
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


void startRFU();
void setoptions(int i, int j, SEXP el, char name[LEN_OPTIONNAME],
		bool isList, bool local);
void getoptions(SEXP sublist, int i, bool local);
void deloptions(bool VARIABLE_IS_NOT_USED local) {
#ifdef DO_PARALLEL
  if (local) RFERROR("'pivot_idx' cannot be freed on a local level");
#endif  
  utilsoption_type *options = WhichOptionList(local);
  FREE(options->solve.pivot_idx);
}

THIS_FILE_ANYSIMD;
EXTERN_SIMD_CHECK(avx2_fctns);
EXTERN_SIMD_CHECK(avx_fctns);
EXTERN_SIMD_CHECK(solve_61);


void loadoptionsRFU() {
  startRFU();
  SetLaMode();

  if (!sseAvail)
    RFERROR("programm does not run on machines that old (not having sse)\n");
  CHECK_THIS_FILE_ANYSIMD;
  CHECK_FILE(avx_fctns);
  CHECK_FILE(avx2_fctns);
  CHECK_FILE(solve_61);
  install_default();

  MEMSET(PIDKEY, 0, PIDMODULUS * sizeof(KEY_type *));
  pid(&parentpid);
  attachRFUoptions((char *) "RandomFieldsUtils",
		   prefixlist, prefixN,
		   allOptions, allOptionsN,
		   setoptions, getoptions,
		   NULL, // final
		   deloptions,
		   NULL, NULL,
		   0, true,
		   GPU_NEEDS, // from configure.ac
		   SIMD_INFO,
		   RFU_VERSION, RFU_VERSION, MEMisALIGNED);

  update_utilsoption();
  KEY_type *KT = KEYT();
  if (KT->global_utils.basic.cores == 0) BUG;

  //finalizeoptions();
}

void  detachoptionsRFU(){
  PIDKEY_DELETE();
  detachRFUoptions(prefixlist, prefixN);
}




#endif // standalone
