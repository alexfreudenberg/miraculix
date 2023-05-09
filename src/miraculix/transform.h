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




#ifndef miraculix_transform_H
#define miraculix_transform_H 1

void haplo2geno(SEXP Code, basic_options *opt);
void haplo2geno(unit_t *M, Long snps, Long individuals, Long lda,
		coding_type oldCoding, int variant, int cores);
coding_type OneSet2G(coding_type coding);
coding_type swapHaploGeno(coding_type coding, bool ToH);
SEXP transpose(SEXP SxI,  option_type *global, utilsoption_type *utils);

SEXP Transform(SEXP SxI, unit_t* SxInt, coding_type codingInfo,
	       int* selSnps, int lenSnps, int *selIndiv, int lenIndiv,
	       option_type *global, utilsoption_type *utils);
#endif


