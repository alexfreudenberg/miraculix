
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


/// sysconf (_SC_NPROCESSORS_ONLN) // number of cores available
// int get_nprocs (void) // dito

#ifndef RFUdef_H
#define RFUdef_H 1

#define STAND_ALONE 1  // for debugging only

//// 1
//#define SCHLATHERS_MACHINE 1
//#define SCHLATHER_DEBUGGING 1


#if ! defined SCHLATHERS_MACHINE && defined SCHLATHER_DEBUGGING
#undef SCHLATHER_DEBUGGING
#else
////  1
#endif

// #define NO_OMP 1

#if ! defined pkg
#define pkg "RandomFieldsUtils" // RFU DELETE
#endif

#endif
