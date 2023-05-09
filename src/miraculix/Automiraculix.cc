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


#include "Automiraculix.h"

const char // muss in separater Zeile sein
// nachfolgende Zeile eingerueckt um 0 oder 2 Zeichen

  *CODING_NAMES[nr_coding_type] = {
  "AutoCoding",
  "1bit",
  "2bit",
  "3bit",
  "plink",
  "5codes8", 
  "orig.plink",  
  "unused7",
  "unused8", "unused9",
  "2bit4",  // unused
  "2bit8",
  "2bit16", // unused
  "2bit32",
  "2bit32-T", 
  "unused15", "unused16",
  "1bit-H", "2bit-H", "1bit8-H", "1bit32-H",
  "1-1-bit32", "2bit64-H", "1-1-bit32-T", "2bit64-HT",
  "dotfile", "filedot", "anyGeno", 
  "5codes8-T", "2bit-T", "orig.plink-T",  "plink-T",
  "unknown"},
  
  *CENTERING_NAMES[nr_centering_type] = {
    "no centering", "row means", "column means", "user"},
  
  *INFO_NAMES[INFO_LAST + 1] = {
    "implementation", "snps", "individuals", "coding", "variant",
    "ldaBitalign", "lda", "genuinehaplo", "missings","expected_size_V",  // 9
    "addr", "align0","memInUnits", "AlignedUnits", "bigendian",   // 14
    "unused", "unused", "unused","unused",
    "current_snps (MoBPS)", // 20
    "unused", "unused", "unused","unused","unused",
    "unused", "unused", "unused","unused", // 29
    "isSNPxInd", "header","DoubledIndividuals", "leadingcolumns", // 33
    "unused","unused", "unused", "unused", "unused", // , 38
    "unused","unused", "unused", "unused", "unused", "unused", // 44
    "RECENTALIGNADDR","BLOCKEDINFO", // 46
    "unused", "unused", "unused","unused", // 50
    "unused","unused", "unused", "unused", "unused", // 55ste
    "unused","unused", "unused", "unused", "unused", "unused", // 61ste
    "Zaehler", "<Last>"},


  *NORMALIZING_NAMES[nr_normalizing_type] = {
    "no normalizing", "according_centering", "global normalization"},
 
  *MIRACULIX_NAMES_OF_NAMES[N_MIRACULIX_NAMES_OF_NAMES] = {
    // add, but NEVER delete or change
    "coding", "centering", "info", "normalizing"
  }

; // 64

