#ifndef MIRACULIX_AUX_H
#define MIRACUKIX_AUX_H 1


#include "local1.h"
#include "error.h"
#include "local3.h"

#include <R.h>
#include <Rinternals.h>



//#include <Rmath.h>
//#include <errno.h>
//#include <R_ext/Complex.h>
//#define MULTIPLEMATCHING -2 //kleinkram.h
//#define NOMATCHING -1
//#define MATCHESINTERNAL -3



extern utilsparam* GLOBAL_UTILS;

typedef char name_type[][MAXCHAR];

void windower_meanC(int *Init, int *Length, int* Step, int *start, int  *ende,
		   double *data, int *Lendata, double *res, int *N);
void windower_minC(int *Init, int *Length, int* Step, int *start, int  *ende,
		  double *data, int *Lendata, double *res, int *N);
void windower_maxC(int *Init, int *Length, int* Step, int *start, int  *ende,
		  double *data, int *Lendata, double *res, int *N);
void windower_medianC(int *Init, int *Length, int* Step, int *start, int *ende,
		     double *data, int *Lendata, double *res, int *N);

void scanC(int *positions, int *length, double *freq, int *minscan,
	  int *maxscan, double *threshold, int *nthres, 
	  int *PER_SNP,
	  int *above_threshold, double *maximum);
void sumscanC(int *positions, int *length, double *freq, int *minscan,
	     int *maxscan, double *threshold, int *nthres, 
	     int *PER_SNP,
	     int *above_threshold, double *maximum);



#endif /* miraculix_aux_h */
