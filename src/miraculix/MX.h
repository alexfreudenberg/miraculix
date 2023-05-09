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


#ifndef miraculix_MX_H
#define miraculix_MX_H 1

#include "MXinfo.h"

#if ! defined MY_VARIANT || ! defined MY_LDABITALIGN
#error MY_VARIANT / MY_LDABITALIGN undefined. Use MXinfo.h instead
#endif


Long fctnLDAalignBitSnps(Long bits, Long LDAbitalign, Long my_LDABITALIGN,
		     coding_type my_CODING);
#define LDAalignBitSnps(bits)						\
  fctnLDAalignBitSnps(bits, LDAbitalign, MY_LDABITALIGN, MY_CODING)


#define LDAalignBlocks(blocks) LDAalignBitSnps((blocks) * BitsPerBlock)


Long static inline Blocks(Long snps) {
  assert(MY_VARIANT <= MY_LDABITALIGN);
  return snps2entities(snps, CodesPerBlock);
}

Long static inline Units(Long snps) { // OK
  return snps2entities(snps, CodesPerUnit); // OK
}
#endif
