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
#include "extern_RFU.h"

 
void startRFU() {
  // utilsoption_type *utils = &OPTIONS;
  union {
    unsigned short a;  
    unsigned char b[2];
  } ab;
  ab.a = 0xFF00;
  bool bigendian = ab.b[0] != 0;
  assert((__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__) xor bigendian);
  OPTIONS.basic.bigendian = bigendian;
  
  volatile uni64 U64;
  U64.u32[bigendian] = 1954;
  U64.u32[!bigendian] = 0x7ff00000;
  OPTIONS.basic.NA = U64.d8[0];
  OPTIONS.basic.NaN = 0.0 / 0.0;

#if defined compatibility_to_C_h  
  ownNA = OPTIONS.basic.NA;
  ownNaN = OPTIONS.basic.NaN;
#endif  

}



