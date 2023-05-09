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


// ACHTUNG: Reihenfolge nicht aendern!
#include "Basic_miraculix.h"
#include "xport_import.h"
#include "MXinfo.h"
#include "intrinsicsCheck.h"
#include "options.h"
#include "mmagpu.h"
#include "utils_miraculix.h"
#include "kleinkram.h"



const char * prefixMlist[prefixMN] =
  {"genetics", "tuning_param", "genetics_messgaes"};

// IMPORTANT: all names of general must have at least 3 letters !!!

const char *genetics_names[geneticsN] = 
  {"digits", "snpcoding", "centered", "normalized", "squared", "haplo_set1",
   "interprete_as_is", "prefer_external_freq"};

const char *tuning_names[tuningN] = 
  {"savehaplo",
   "oldversion", // for debugging only !
   "variant", //variant
   "logSnpGroupSize", // min 12
   "indivGroupSize", // min 1
   "smoothSnpGroupSize"
   "meanVsubstract", "meanSxIsubstract"
   "floatPrecision",
   "missingsFully0",
   "miniRows", "miniCols", "use_gpu_for_vector",
   "add_tranposed"
  };


const char *messages[messagesN] = 
  {"warn_address"}; 


option_type OPTIONS_MIRACULIX = { // OK
  genetics_START,
  tuning_START,
  messages_START
};


const char **allMoptions[prefixMN] = {genetics_names,tuning_names,messages};
int allMoptionsN[prefixMN] = {geneticsN,tuningN,messagesN};


int main_variant(int variant) {
  return variant < 128 ? variant: variant - variant % 128;
  //  return variant < 128 ? variant
  //  : variant - variant % (variant < VARIANT_GPU ? 256 : 128);
}
 
bool exists_variant(coding_type coding, int variant, bool indeed,
		    bool efficient) {
   assert(FirstUserHaploCoding == 14 && LastUserHaploCoding == 18);
 //  printf("var = %d %d %d\n", variant, coding, indeed);
  if (variant == VARIANT_DEFAULT) return true;
  if (variant < 0 || variant > 32 *(NvariantsSHR5 - 1) || (variant % 32 != 0)
      || !variantsSHR5[variant >> 5])
    return false;
  
  bool ans;
  switch(coding) {
  case OneBitGeno : ans = variant < VARIANT_GPU; break;
  case TwoBitGeno :
    ans = variant <= 512 + VARIANT_A &&
      (variant <= 64 || (variant % 128) <= VARIANT_A)
#if defined USEGPU
      || variant == VARIANT_GPU 
#endif
      ; break;
  case ThreeBit : ans = variant == 64; break;
  case FourBit : ans = false; break;
  case OneByteGeno : ans = variant == 32;  break;
  case TwoByte : ans = false; break;
  case FourByteGeno : ans = variant == 32 || variant == VARIANT_R ||
      variant == VARIANT_R + VARIANT_A
#if defined USEGPU
      || variant == VARIANT_GPU 
#endif
      ; break;
  case FourByte_Tr : ans = false; break;
  case OneBitHaplo : ans = false; break;
  case TwoBitHaplo: ans = variant == 32; break;
  case OneByteHaplo: ans = variant == 32; break;
  case FourByteHaplo: ans = variant == 32; break;
  case FourByteSingleBit: ans = variant == 32; break;
  case EightByteHaplo: ans = variant == 32; break;
  case FourByteSingleBit_Tr: ans = variant == 32; break;
  case EightByteHaplo_Tr: ans = variant == 32; break;
  default : BUG;
  }

  //  printf("ok=%d (%d) %d %d\n", ans, efficient, coding, variant);

  // printf("%s %d indeed=%d: ans=%d -> %d; #!indeed=%d ok=%d\n",	 CODING_NAMES[coding], variant, indeed, ans,	 check_variant(coding, variant, efficient, false),  !ans || !indeed,	 (int) check_variant(coding, variant, efficient, false) >= (int) VARIANT_DEFAULT);
  
  if (!ans || !indeed) return ans;
  
  return check_variant(coding, variant, efficient, false) >= VARIANT_DEFAULT;
}


bool exists_tiling(coding_type coding, int variant, bool tiling) {
  // function assumes that combination coding/variant is valid!
 assert(FirstUserHaploCoding == 14 && LastUserHaploCoding == 18);
  switch(coding) {
  case OneBitGeno : return true;
  case TwoBitGeno :
    return (variant >= 128 && variant < VARIANT_GPU) ||
      ((variant < 128 || variant == VARIANT_GPU) && !tiling);
  case ThreeBit : return !tiling;
  case FourBit : return false; 
  case OneByteGeno : return false; 
  case TwoByte : return false; 
  case FourByteGeno : return !tiling;
  case FourByte_Tr : return false; 
  case OneBitHaplo : return false; 
  case TwoBitHaplo: return false; 
  case OneByteHaplo: return false; 
  case FourByteHaplo: return false; 
  case FourByteSingleBit: return false; 
  case EightByteHaplo: return false; 
  case FourByteSingleBit_Tr: return false; 
  case EightByteHaplo_Tr: return false; 
  default : BUG;
  }
  BUG;
  return false;  
 }
  


bool exists_crossprod(coding_type m) {
   assert(FourByteGeno == 13);
   return (m >= FirstGenoCoding && m <= ThreeBit) || m == FourByteGeno;
}

bool exists_allelefreq(coding_type m) {
  assert(FiveCodes == 5);
  return  (m >= TwoBitGeno && m <= FiveCodes) ||
    m == FiveCodesTransposed || m == TwoBitGenoTransposed ||
    m == OneByteGeno || m == FourByteGeno ||
    (m >= OneByteHaplo && m <= EightByteHaplo_Tr);
}


bool exists_internal_coding(coding_type m) {
  assert(FourByteGeno == 13);
  return (m >= OneBitGeno && m <= FourByteGeno &&
	  m != FourBit && m != TwoByte) ||
    (m >= OneBitHaplo && m <= OneByteHaplo);
}


bool exists_user_coding(coding_type m) {
  assert(EightByteHaplo_Tr == 24);
  return m == OneByteGeno || m == FourByteGeno ||
    (m >= OneByteHaplo && m <= EightByteHaplo_Tr);
}

bool exists_coding(coding_type m) {
  return exists_internal_coding(m) || exists_user_coding(m);
}

