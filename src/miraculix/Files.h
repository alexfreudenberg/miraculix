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



// this is mainly for plinkbinary files



void CheckSnpsLessIndiv(Long snps, Long individuals, bool isSNPxInd,
			basic_options *opt);

typedef void (*TrafoFileValues_t) (unit_t *CodeFile, Long blocks);

typedef void (*coding_sexp_t)(unit_t *, Long, Long, Long, Long, Long,
			      basic_options *,double *, SEXP);

SEXP file_intern(SEXP file, coding_type coding, int variant,
		 option_type *global, utilsoption_type *utils,
		 bool only_complete,
		 SEXP G);

SEXP file_binary(char *file, char *file_coding, int header, bool isSNPxInd, 
		 // output
		 Long snps,
		 Long individuals,
		 Long codesperblock, // OK
		 coding_type coding, int variant,
		 Long  LDAbitalign,
		 coding_sexp_t coding_sexp,
		 option_type *global, utilsoption_type *utils,
		 SEXP G
		 );



/*
  R = AND( N1E1, a);					\
  L = AND( E1N1, a); 		\
  tmp = XOR( E1N1, L); 				\
  tmp = SHR64( tmp, 1); 				\
  tmp = AND( tmp, R);		\
  if (ANY(tmp)) ERR0("missing value detected.");			\
  tmp = SHR64( L, 1); 				\
  tmp = XOR( L, R);						\
  R = SHL64( R, 1); 				\
  L = AND( L, R);							\
  tmp = OR( tmp, L);   						\
  STORE(block + j,  tmp);						\
*/

/* Plink:
   00 -> 00
   01 : missing
   10 -> 01
   11 -> 10
*/
#define FileBinary(NR)						\
  void TrafoFileValues_##NR(unit_t *CodeFile, Long blocks) {	\
    BlockType0* code = (BlockType0*) CodeFile;			\
    for (Long j=0; j<blocks; j++) {				\
      const BlockType a = LOAD(code + j);			\
      const BlockType						\
	R = AND( N1E1, a),					\
	L = AND( E1N1, a), /* 1* -> 10, 0*->00 */		\
	c1 = XOR( E1N1, L), /* 00, 10 */			\
	c2 = SHR64( c1, 1), /* 00, 01 */			\
	c3 = AND( c2, R); /* was the code 01 anywhere ? */	\
      if (ANY(c3)) ERR0("missing value detected.");		\
      const BlockType						\
	d1 = SHR64( L, 1), /* right digit */			\
	d2 = XOR(d1, R),					\
	R1 = SHL64( R, 1), /* left digit */			\
	L1 = AND( L, R1),					\
	d3 = OR( d2, L1);					\
      STORE(code + j,  d3);					\
    }								\
  }
     
		

#define UPDATE	/*//printf("OXXAAK\n");*/				\
  assert(coding_sexp != NULL);						\
  if (isSNPxInd) { /*//printf("OK %d %d\n", ( r + 1) % nrow_matrix, coding_sexp==NULL);*/ \
    if ((r + 1) % nrow_matrix == 0) {					\
      /* always delivering matrix as snps x individuals */		\
      coding_sexp (matrix, matrixLDA,					\
		   r + 1-nrow_matrix, r + 1, 0, individuals, 		\
		   opt, dG, Ans);					\
      rowidx = -plusrow;						\
      MEMSET(matrix, 0, matrix_size * BytesPerUnit); /* // OK */	\
    }									\
  } else { /*// printf("OXXK\n");*/					\
    coding_sexp(matrix, matrixLDA,					\
		0, snps, r, r + 1,					\
		opt, dG,Ans);					\
    rowidx = -plusrow;							\
    MEMSET(matrix, 0, matrix_size * BytesPerUnit); /* // OK */			\
  }


#define CLOSE	/* read from file rowwise */				\
  if (isSNPxInd) {							\
    Long remainder = r % nrow_matrix;					\
    if (remainder != 0) {						\
      coding_sexp(matrix, matrixLDA,					\
		  (Long)(r - remainder), r ,0, individuals, 		\
		  opt, dG,  Ans);				\
    }									\
  } /* KEIN else  */							\
  fclose(fp);								\
  FREE(matrix); 
  
