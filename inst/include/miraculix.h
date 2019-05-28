/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2018 -- 2018  Martin Schlather

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


#ifndef miraculix_main_H
#define miraculix_main_H 1

//
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif 
  //  void testint(int *a, int*b, int*N, int*add);
  //  void testdouble(double *a, double*b, int*N, int*add);

 // general
  SEXP RFoptions(SEXP options);
  void RelaxUnknownRFoption(int *relax);

  void attachmiraculix();
  void detachmiraculix();
 
  // scan
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


  // windower
  SEXP windower_mean(SEXP Init, SEXP Length, SEXP Step, SEXP start, SEXP ende,
		     SEXP data, SEXP Lendata, SEXP N);
  SEXP windower_min(SEXP Init, SEXP Length, SEXP Step, SEXP start, SEXP ende,
		    SEXP data, SEXP Lendata, SEXP N);
  SEXP windower_max(SEXP Init, SEXP Length, SEXP Step, SEXP start, SEXP ende,
		    SEXP data, SEXP Lendata, SEXP N);
  SEXP windower_median(SEXP Init, SEXP Length, SEXP Step, SEXP start, SEXP ende,
		       SEXP data, SEXP Lendata, SEXP N);

  // Verwandtschaft
  SEXP vector012matrix(SEXP vector, SEXP matrix);
  SEXP matrixvector012(SEXP vector, SEXP matrix);

  
  // Hamming for SSE2 and SSSE3
  SEXP sse_check();
  SEXP matrixH_get(SEXP M);  
  SEXP matrix_codingH(SEXP MT2);
  SEXP matrix_mult(SEXP M);

  
  // large table coding with 2 or 3 bits
  SEXP matrix_coding(SEXP M) ;  //
  SEXP file_coding(SEXP file);
  SEXP MmPsqSNPs(SEXP SNPxIndiv, SEXP compressed_start, SEXP compressed_end);
  SEXP MmPsqIndividualsParallel(SEXP z);
  SEXP getRelmatrixIndividuals(SEXP z, SEXP SNPxIndiv);
  SEXP getRelmatrixSNP(SEXP z, SEXP SNPxIndiv);
  SEXP file_get(SEXP file);
  SEXP file_dot(SEXP file, SEXP g);
  SEXP dot_file(SEXP file, SEXP g);
  SEXP matrix_get(SEXP M);  
  
  SEXP codeOrigins(SEXP M);
  SEXP decodeOrigins(SEXP M, SEXP line);

  // shot table coding
  SEXP fillSNPmatrix(SEXP Z, SEXP Idx, SEXP V) ;
  SEXP createSNPmatrix(SEXP SNPS, SEXP Individuals) ;  
  SEXP relationshipMatrix(SEXP CGM) ;
  SEXP decodeGeno(SEXP CM);
  SEXP decodeNthGeno(SEXP CM, SEXP NN);
  SEXP decodeHaplo(SEXP CM) ;
  SEXP haplo2geno(SEXP Haplo);
  SEXP codeGeno(SEXP G) ;
  SEXP zeroNthGeno(SEXP CM, SEXP NN);
  SEXP copyGeno(SEXP CM);
  SEXP codeHaplo(SEXP M) ;// 2 x Haplo - Matrix, as vector gespeichert
  SEXP codeHaplo2(SEXP M1, SEXP M2);
  SEXP vectorGeno(SEXP V, SEXP Z);
  SEXP genoVector(SEXP Z, SEXP V);
  SEXP unlock(SEXP X);
  SEXP dolocking(SEXP Do);
  SEXP computeSNPS(SEXP population, SEXP Generation, SEXP Sex, SEXP Nr,
		   SEXP From_SNP, SEXP To_SNP, SEXP Select, SEXP Geno);
  SEXP allele_freq(SEXP CGM);
  SEXP Debug();
  SEXP StopDebug();
  SEXP compute(SEXP population, SEXP Generation, SEXP Sex, SEXP Nr,
	       SEXP tau, SEXP vec, SEXP b, SEXP Select, SEXP matrix_return);
  SEXP solveRelMat(SEXP R, SEXP tau, SEXP vec, SEXP b, SEXP destroy);
  SEXP substract_centered(SEXP SnpXindiv);
  SEXP get_centered();
  
  SEXP failure(SEXP file);

  SEXP init_parallelcompute();

#ifdef __cplusplus
}
#endif

#endif /* RF_simu_PUBLIC_H*/
