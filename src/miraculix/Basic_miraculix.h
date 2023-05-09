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

#ifndef basic_miraculix_H
#define basic_miraculix_H 1

#include "def.h"
#include "intrinsics.h"
#include "compatibility.general.h"
#include "compatibility.SEXP.h"

#include "Basic_RandomFieldsUtils.h"

//typedef Uint unit_t;

/*
typedef unsigned int unit_t;
typedef unsigned int _4Byte;
typedef unsigned char _1Byte;
*/

/*
typedef short unsigned int unit_t;
typedef Uint _4Byte;
typedef Ulong  _1Byte;
*/

/*
typedef short unsigned int unit_t;
typedef Uint _1Byte;
typedef Ulong  _4Byte;
*/

/*
typedef Uint unit_t;
typedef short unsigned int _1Byte;
typedef Ulong  _4Byte;
*/

/*
typedef int unit_t;
typedef Ulong _1Byte;
typedef short unsigned int  _4Byte;

*/

/*
typedef Long unit_t;
typedef int _1Byte;
typedef short unsigned int  _4Byte;
*/

/*
typedef Long unit_t;
typedef short unsigned int _1Byte;
typedef  int  _4Byte;
*/

typedef unsigned int unit_t; // general code
typedef unsigned char _1Byte;         // half user OneByte
typedef unsigned int  _4Byte;// user 4-Byte

#endif
