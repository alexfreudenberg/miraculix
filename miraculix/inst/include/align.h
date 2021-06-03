
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

#ifndef miraculix_align_H
#define miraculix_align_H 1

#define FROMINPUT  *(pX++)
#define FROMHAPLO GetHaplo(pX, s)

#if defined CodesPerUnit
Ulong static inline Blocks(Ulong snps) {				
  return 1L + (snps - 1L) / CodesPerBlock; }				
Ulong static inline Units(Ulong snps) {				
  return 1L + (snps - 1L) / CodesPerUnit; }			       
#define Align(CM, nr) AlignBase(CM, nr, BytesPerBlock, true)
#define AlignTest(CM, nr, test)	AlignBase(CM, nr, BytesPerBlock, test)
#endif 

#define Align256(CM, nr) AlignBase(CM, nr, BytesPerBlock256, true)
#define Align256test(CM, nr, test) AlignBase(CM, nr, BytesPerBlock256, test)

#endif
