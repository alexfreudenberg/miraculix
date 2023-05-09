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

#ifndef RandomFieldsUtilsxport_H
#define RandomFieldsUtilsxport_H 1

typedef
struct KEY_type KEY_type;
struct KEY_type {
  KEY_type *next;
  utilsoption_type global_utils;
  int pid,  visitingpid;
  bool ok, doshow;
  errorstring_type error_location;

  int *ToIntDummy;
  int ToIntN, ToRealN ;
  double *ToRealDummy;

  whittle_work_type whittle;
};
extern KEY_type *PIDKEY[PIDMODULUS];
KEY_type *KEYT();

typedef
struct option_type option_type;
utilsoption_type *WhichOptionList(bool local);
void PIDKEY_DELETE();

extern const char *R_TYPE_NAMES[LAST_R_TYPE_NAME + 1];

#endif
