/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de

 Copyright (C) 2017 -- 2019 Martin Schlather

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

#ifndef AutoMiraculix_H
#define AutoMiraculix_H 1
#include <R.h>

typedef enum snpcoding {AutoCoding,  // 0
			NoSNPcodingR, // none needed, 1
			ThreeBit, // none needed, 2
			Hamming2, // SSE2 3
			Hamming3, // SSSE3, 4
			NoSNPcoding, // AVX2 needed, 5
			TwoBit,   // none needed 6
			Multiply, // 7
			Packed,  // SSE2 8
			Shuffle,  // SSSE3 9
			Multiply256, //10
			Packed256, // 11
			Shuffle256, // or AVX2 12
			MMAGPU,
			CaseCount, // ab AVX512
			unused15,
			unused16, 
			unused17, 
			unused18, 
			unused19,
			unused20,
			unused21,
			unused22,
			unused23,
			unused24,
			unused25,
			unused26,
			unused27,
			unused28,
			unused29,
			Haplo, // 30; see MX.h! before changing!!
			UnknownSNPcoding // 31
} snpcoding;

#define nr_snpcoding (UnknownSNPcoding + 1)
#define FirstMoBPSmethod TwoBit
#define LastMoBPSmethod MMAGPU
#define FirstGenuineMethod ThreeBit
#define LastGenuineMethod LastMoBPSmethod


// WHAT just doubles the information that is also available through
// the class information, but easier to access on the C level for
// historical reasons. To be deleted maybe somewhen.


// Coding of the attribute "information"
// !!!! ACHTUNG ZAHLEN MUESSEN DIE GLEICHEN BLEIBEN !!!!

#define CURRENT_VERSION 3

#define VERSION 0
#define SNPS 1 // Wert darf auf keinen Fall geaendert werden
#define INDIVIDUALS 2 // Wert darf auf keinen Fall geaendert werden
#define ADDR0 3
#define ADDR1 4 // Achtung! Zweiter Teil von ADDR !!
#define ALIGNADDR0 5
#define ALIGNADDR1 6
#define SUMGENO 7
#define SUMGENO_E9 8
#define METHOD 9
#define ALIGNMENT 10
#define SNPxIND 11 // INFO_INDIV_PER_COL, i.e. 'percolumn'
#define BITSPERCODE 12
#define BYTESPERBLOCK 13
#define CODESPERBLOCK 14
#define HEADER 15
#define DOUBLEINDIV 16
#define LEADINGCOL 17
#define MEMinUNITS0 18  // total genuine memory used to store all matrices --
// neither memory for allignment is not included, nor memory for info
#define MEMinUNITS1 19  
#define ALIGNEDUNITS0 20  // memory needed including alignment
#define ALIGNEDUNITS1 21
// ACHTUNG: gegebenenfalls haplogeno.cc:copyGeno aendern!!
#define UNITSPERINDIV 22
#define INFO_GENUINELY_LAST UNITSPERINDIV
#define RECENTALIGNADDR0 23
#define RECENTALIGNADDR1 24


#define BLOCKEDINFO 61
#define ZAEHLER 62
#define INFO_LAST 63

#define CURRENT_SNPS HEADER // used in MOBPS, which doesn't use HEADER

#define GENOMICMATRIX "genomicmatrix"
#define HAPLOMATRIX "haplomatrix"
#define ORIGINVECTOR "origindata"

extern const char *SNPCODING_NAMES[nr_snpcoding],
  *INFO_NAMES[INFO_LAST + 1];
#endif
