
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2018 -- 2019  Martin Schlather

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/

#include "AutoMiraculix.h"

const char // muss in separater Zeile sein
// nachfolgende Zeile eingerueckt um 0 oder 2 Zeichen
*SNPCODING_NAMES[nr_snpcoding] = {
  "AutoCoding", "NoSNPcodingR",
  "ThreeBit",
  "Hamming2",
   "Hamming3",
  "NoSNPcoding",
  "TwoBit",
  "Multiply","Packed",
  "Shuffle",
  "Multiply256","Packed256", "Shuffle256",
  "MMAGPU", "CaseCount",  
  "unused", "unused", "unused",
  "unused", "unused", "unused", "unused", "unused",
  "unused", "unused", "unused", "unused", "unused",
  "Haplo", "unknown"},
  
  *INFO_NAMES[INFO_LAST + 1] = {
  "version", "snps", "individuals", "addr0", "addr1",
  "align0", "align1", "sumgeno", "sumgenoE9", "method",
  "alignment", "isSNPxInd", "bitspercode", "bytesperblock", "codesperblock",
  "header", "DoubledIndividuals", "leadingcolumns", "memInUnits0","meminUnits1",
  "AlignedUnits0", "AlignedUnits1", "unitsperindiv", "unused", "unused",
  "unused", "unused", "unused","unused", "unused",  // 30
  "unused","unused", "unused", "unused", "unused",
  "unused","unused", "unused", "unused", "unused",  // 40
  "unused","unused", "unused", "unused", "unused",
  "unused","unused", "unused", "unused", "unused", // 50
  "unused","unused", "unused", "unused", "unused",
  "unused","unused", "unused", "unused", "unused", // 60
  "unused","unused", "unused", "unused"}; // 64
