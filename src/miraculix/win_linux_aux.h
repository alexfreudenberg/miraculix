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



#ifndef WIN_LINUX_AUX_H
#define WIN_LINUX_AUX_H 1

uint32_t cpuid_info(int Blatt, int Register);//MINGWCPUID, WINCPUID, LINUXCPUID

#ifdef __cplusplus
extern "C" {
#endif
  void sleepMilli(int *milli);
  void sleepMicro(int *milli);
  void pid(int *i);
  void hostname(char **h, int *i);
  bool
    parallel();
#ifdef __cplusplus
}
#endif


#endif /* WIN_LINUX_AUX_H */


