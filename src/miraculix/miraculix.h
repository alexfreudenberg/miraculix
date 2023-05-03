
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

#ifndef miraculix_main_H
#define miraculix_main_H 1

#ifdef __cplusplus
extern "C" {
#endif 

 // general
  //  SEXP RFoptions(SEXP options);
  void loadoptions(int *n);
  void detachoptions();
  SEXP copyoptions();

  SEXP existsVariant(SEXP Coding, SEXP Variant, SEXP Indeed);
  SEXP existsTiling(SEXP Coding, SEXP Variant, SEXP Tiling);
  SEXP existsCrossprod(SEXP Coding);
  SEXP existsAllelefreq(SEXP Coding);
  SEXP existsCoding(SEXP Coding, SEXP Internal);
  SEXP get_centered(); 
 
  // scan / windower
  SEXP scan(SEXP positions, SEXP length, SEXP freq, 
	    SEXP minscan,  SEXP maxscan, SEXP threshold, SEXP nthres,
	    SEXP PER_SNP,
	    SEXP above_threshold, SEXP maximum);
  
  SEXP sumscan(SEXP positions, SEXP length, SEXP freq, 
	       SEXP minscan,  SEXP maxscan, SEXP threshold, SEXP nthres,
	       SEXP PER_SNP,
	       SEXP above_threshold, SEXP maximum);

  SEXP collect_scan(SEXP positions, SEXP length, SEXP freq, 
		    SEXP minscan,  SEXP maxscan, SEXP threshold, SEXP nthres,
		    SEXP PER_SNP, 
		    // additional return values:
		    SEXP areas,  SEXP value
		    );
  SEXP collect_scan2(SEXP positions, SEXP length, SEXP freq, 
		     SEXP minscan,  SEXP maxscan, SEXP threshold, SEXP nthres,
		     SEXP PER_SNP, SEXP max_intervals,
		     SEXP max_max_basepair_dist, SEXP exclude_negative,
		     // additional return values:
		     SEXP above_threshold, SEXP maximum		     
		     );

  SEXP windower(SEXP Init, SEXP, SEXP Length, SEXP Step, SEXP start, SEXP ende,
		     SEXP data, SEXP Lendata, SEXP N);


  // creation of compressed SNP matrix  
  SEXP createSNPmatrix(SEXP SNPs, SEXP Individuals) ;  
  SEXP fillSNPmatrix(SEXP Z, SEXP Idx, SEXP V) ;
  SEXP zeroGeno(SEXP Code, SEXP SNPs, SEXP Indiv, SEXP copy);
  SEXP rhaplomatrix(SEXP freq, SEXP freq2,  SEXP individuals, SEXP coding);
  SEXP Transform(SEXP M, SEXP Coding, SEXP SelSNPs, SEXP SelIndiv);
  

  // vector-matrix operation with mixed SNP non-SNP arguments
  SEXP vector012matrix(SEXP vector, SEXP matrix);
  SEXP matrixvector012(SEXP vector, SEXP matrix);
  SEXP vectorGeno(SEXP V, SEXP Z);
  SEXP genoVector(SEXP Z, SEXP V);

 
  // mathematical operations on SNP matrix
  SEXP transpose(SEXP x);
  SEXP crossprod(SEXP M);
  SEXP allele_freq(SEXP SxI);
  SEXP substract_centered(SEXP SnpXindiv);


  // mathematical operations including relationship matrix
  SEXP solveRelMat(SEXP R, SEXP tau, SEXP vec, SEXP b, SEXP destroy);

  
  // MoBPS
  SEXP codeOrigins(SEXP M);
  SEXP decodeOrigins(SEXP M, SEXP line);
  SEXP computeSNPS(SEXP population, SEXP Generation, SEXP Sex, SEXP Nr,
		   SEXP From_SNP, SEXP To_SNP, SEXP Select, SEXP Geno);
  SEXP compute(SEXP population, SEXP Generation, SEXP Sex, SEXP Nr,
	       SEXP tau, SEXP vec, SEXP b, SEXP Select, SEXP matrix_return);
  SEXP MiraculixOptions(SEXP options);
  SEXP Debug();
  SEXP StopDebug();

  // others
  SEXP crossprodInt(SEXP X, SEXP Y, SEXP mode); // necessary? joint with 4byte?
 

#ifdef __cplusplus
}
#endif

#endif /* RF_simu_PUBLIC_H*/
