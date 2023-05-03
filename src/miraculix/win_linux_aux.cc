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



#include "def_rfu.h"
#include "intrinsics.h"

#if defined MSDOS_WINDOWS
#define VC_EXTRALEAN
#include <windows.h>
#endif


// achtung! windows.h zusammen mit <Rmath.h oder R.graphics>
// gibt warnung, da ERROR mehrfach definiert !
// deshalb auch in auxiliary.h nicht basic.h einbinden // obsolette ?!!
#include <unistd.h>
#if defined compatibility_to_R_h
#include <Rinternals.h>
#endif
#include "win_linux_aux.h"


 
void sleepMilli(int *milli) {
#if defined MSDOS_WINDOWS
  Sleep((long) *milli); // along
#else 
  usleep((useconds_t) (1000 * (unsigned long) *milli));// along
#endif
}

void sleepMicro(int *micro) {
#if defined MSDOS_WINDOWS
  Sleep((long) ((*micro + 500) / 1000));// along
#else
  usleep((useconds_t) *micro);
#endif
}

void pid(int *i)  {
#if defined MSDOS_WINDOWS
  *i = _getpid();
#else
  *i = getpid(); 
#endif
}

int parentpid=0;
bool
  parallel() {
    int mypid;
    pid(&mypid);
    return mypid != parentpid;
  }


void hostname(char **h, int *i){
#if defined MSDOS_WINDOWS
  *h[0]=0;
#else
  gethostname(*h, *i);
#endif
}  

#define VARIABLE_IS_NOT_USED __attribute__ ((unused))
uint32_t cpuid_info(int VARIABLE_IS_NOT_USED Blatt,
		    int VARIABLE_IS_NOT_USED Register) {
#if defined MINGWCPUID
   uint32_t s[4];						     
   __cpuid(Blatt, s[0], s[1], s[2], s[3]);			    
   return s[Register];
#elif defined WINCPUID
  uint32_t s[4];							
  __cpuid((int *)s, (int) Blatt);
  return s[Register];
#elif defined LINUXCPUID
  uint32_t s[4];							
  asm volatile							
    ("cpuid": "=a"(s[0]), "=b"(s[1]),"=c"(s[2]),			
     "=d"(s[3]):"a"(Blatt),"c"(0));					
  return s[Register]; 
#else
   return 0;
#endif
}
