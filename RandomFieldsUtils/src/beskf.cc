/*gcc beskf.c -o [name] -lm -mavx*/

/*Modified Bessel Function*/

/*From R-3.6.1/src/nmath/b essel_k.c*/

//#include <t ime.h>

#include <stdio.h>
//#include <stdlib.h>
//#include <m ath.h>
#include <float.h> /*DBL_MIN*/
#include <immintrin.h> /*AVX -mavx*/

#include <Rmath.h>
#include "intrinsics.h"
#include "General_utils.h"

#define xmax_BESS_K	705.342 /*From bessel.h*/
#define sqxmin_BESS_K	1.49e-154 /*From bessel.h*/

#define M_SQRT_2dPI	0.797884560802865355879892119869 /*From Rmath.h*/

#define min0(x, y) (((x) <= (y)) ? (x) : (y))
#define max0(x, y) (((x) <= (y)) ? (y) : (x))


/*void bes_k_simd (double *xv, double alpha, int sx, double *yv);

int main (void)
{
	int i, len = 1000000;
	double alpha = 1.1;
	double *x= m alloc(len *sizeof(double)),
			*yv = m alloc(len *sizeof(double)),
			timetaken;
	clock_t start, end;
	
	for (i = 0; i < len; i++){
		x[i] = (double)rand()/RAND_MAX*5.0;
		yv[i] = 0.;
	}
	start = clock();
	bes_k_simd(x,alpha,len,yv);
	end = clock();
	timetaken = (double) (end - start);
//	printf("Time: %f\n",timetaken);
	f ree(x);
	f ree(yv);
	return 0;
}*/

void bes_k_simd (double VARIABLE_IS_NOT_USED *xv,
		 double VARIABLE_IS_NOT_USED Nu,
		 int VARIABLE_IS_NOT_USED  sx,
		 double VARIABLE_IS_NOT_USED *yv)
{
  ERR("not yet programmed");
}
