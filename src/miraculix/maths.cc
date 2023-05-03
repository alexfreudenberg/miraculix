
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



#include "Basic_RFUlocal.h"
#include "compatibility.lapack.h"
#include "zzz_RFU.h"
#include "utils.h"

#if defined compatibility_to_R_h
#include "xport_import_RFU.h"

#else

double bessel_k_ex(double VARIABLE_IS_NOT_USED y,
		   double VARIABLE_IS_NOT_USED nu,
		   double VARIABLE_IS_NOT_USED logscale,
		   double VARIABLE_IS_NOT_USED *bk) {
  ERR0("besselK not available in standalone version");
  return RF_NAN;
}
double pgamma(double x, double shape, double scale, int lower, int logarithm){
/*
f(x, shape, scale, lower, log) {
  x <- x /  scale;
  f = 1.0 + shape;
  quotient = 1.0;
  sum = 0.0;
  res =  shape * log(x) - x - lgamma(f);
  while (quotient >= 1.0) {
      quotient = quotient * x / f;
      sum = sum + quotient;
      f = f + 1.0;
  }
  while (quotient >= eps * sum) {
    quotient = quaotion * x / f;
    sum = sum +  quotient;
    f = f +  1.0;
  }
  res += LOG1P(sum);
  if (lower)
  if (log 
}

*/
  double eps = 1e-17;
  // double logthreshold = -690;
  x /= scale;
  //  double upper_logupperbound = 1.0 - shape - (x + 1) / (shape + x) * x
    //    + shape * LOG(shape + x) - LOG1P (x);
  //  if (lower) {
  //    if (logarithm && upper_logthreshold < -40) return 1.0;
  //  } if (upper_logupperbound < logthreshold) 
  //    if (logarithm && upper_logthreshold < -40) return 1.0;
  //  }
  double
    f = 1.0 + shape,
    quotient = 1.0,
    sum = 0.0;
  double res =  shape * LOG(x) - x - lgamma(f);
  while (quotient >= 1.0) {
      quotient *= x / f;
      sum += quotient;
      f += 1.0;
  }
  while (quotient >= eps * sum) {
    quotient *= x / f;
    sum += quotient;
    f += 1.0;
  }
  res += LOG1P(sum);
  if (lower && logarithm) return res;
  res = EXP(res);
  if (lower) return res;
  res = 1.- res;
  if (logarithm) return log(res);
  return res;
}

#endif

double struve_intern(double x, double nu, double factor_Sign, bool expscaled)
{ 
 if ((x == 0.0) && (nu>-1.0)) return 0.0;
  if (x <= 0.0) return RF_NA; // not programmed yet
  double exp_dummy,
     dummy = 0.0, 
     logx = 2.0 * LOG(0.5 * x), 
     x1 = 1.5, 
     x2 = nu + 1.5,
     value = 1.0, 
     fsign = factor_Sign,
    epsilon=1e-20;

  
   do {
     dummy += logx - LOG(x1) - LOG(FABS(x2));
     exp_dummy = EXP(dummy);
     value += (1 - 2 * (x2 < 0))  * fsign * exp_dummy;
     //  printf("%10g %10g %10g %10g\n", value, fsign, x1, x2);
     x1 += 1.0;
     x2 += 1.0;
     fsign = factor_Sign * fsign; 
   } while (exp_dummy > FABS(value) * epsilon);

   x1 = 1.5;
   x2 = nu + 1.5;
   if (x2 > 0.0) { 
     dummy = (nu + 1.0) * 0.5 * logx - LOGGAMMAFN(x1) - LOGGAMMAFN(x2);
     if (expscaled) dummy -= x;
     value *= EXP(dummy);
   } else {
     //if ( (double) ((int) (x1-0.5)) != x1-0.5 ) return RF_NA;
     value *= POW(0.5 * x, nu + 1.0) / (GAMMAFN(x1) * GAMMAFN(x2));
     if (expscaled) value *= EXP(-x);
   }

  return value;
}

//void Struv eH(double *x, double *nu) {*x=struv e(*x, *nu, -1.0, false);}
//void Struv eL(double *x, double *nu, int * expScaled) {
//  *x=struv e(*x, *nu, 1.0, (bool) *expScaled);
//}
double StruveH(double x, double nu) {return struve_intern(x, nu, -1.0, false);}
double StruveL(double x, double nu, bool expScaled) {
 return struve_intern(x, nu, 1.0, expScaled);
}





//void I0ML0(double *x, int *n) {
//  int i;
//  for (i=0; i<*n; i++) x[i] = I0mL0(x[i]);
//} 


double I0mL0(double x){
  /* Bessel I_0 - Struve L_0 for non-negative arguments x */
  /* see J. MacLeod, Chebyshev expansions for modified {S}truve and 
     related functions, Mathematics of Computation, 60, 735-747, 1993 */
  static double g2[24] = {
	0.52468736791485599138e0,
	-0.35612460699650586196e0,
	0.20487202864009927687e0,
	-0.10418640520402693629e0,
	0.4634211095548429228e-1,
	-0.1790587192403498630e-1,
	0.597968695481143177e-2,
	-0.171777547693565429e-2,
	0.42204654469171422e-3,
	-0.8796178522094125e-4,
	0.1535434234869223e-4,
	-0.219780769584743e-5,
	0.24820683936666e-6,
	-0.2032706035607e-7,
	0.90984198421e-9,
	0.2561793929e-10,
	-0.710609790e-11,
	0.32716960e-12,
	0.2300215e-13,
	-0.292109e-14,
	-0.3566e-16,
	0.1832e-16,
	-0.10e-18,
	-0.11e-18
  };
  static double g3[24] = {
	2.00326510241160643125e0,
	0.195206851576492081e-2,
	0.38239523569908328e-3,
	0.7534280817054436e-4,
	0.1495957655897078e-4,
	0.299940531210557e-5,
	0.60769604822459e-6,
	0.12399495544506e-6,
	0.2523262552649e-7,
	0.504634857332e-8,
	0.97913236230e-9,
	0.18389115241e-9,
	0.3376309278e-10,
	0.611179703e-11,
	0.108472972e-11,
	0.18861271e-12,
	0.3280345e-13,
	0.565647e-14,
	0.93300e-15,
	0.15881e-15,
	0.2791e-16,
	0.389e-17,
	0.70e-18,
	0.16e-18
  };
  double r, x2, ac;
  int i;
  
  if (x < 0.0) {return RF_NA;}
  if (x < 16.0) {
    r = 0.5 * g2[0];
    ac = ACOS((6.0 * x - 40.0) / (x + 40.0));
    for (i=1; i<24; i++) {
      r += g2[i] * COS(i * ac);
    }
  } else {
    r = 0.5 * g3[0];
    x2 = x * x;
    ac = ACOS((800.0 - x2) / (288.0 + x2));
    for (i=1; i<24; i++) {
      r += g3[i] * COS(i * ac);
    }
    r *= T_PI /* 2/pi */ / x;
  }
  return r;
}




/* Gausian model */
double Gauss(double x) {
  return EXP(- x * x);
  //  printf("%10g %10g\n", *x, *v);
}
double logGauss(double x) {
  return - x * x;
}
double DGauss(double y) {
  return -2.0 * y * EXP(- y * y);
}
double DDGauss(double x) {
  double y = x * x; 
  return (4.0 * y - 2.0)* EXP(- y);
}
double D3Gauss(double x) {
  double y = x * x; 
  return x * (12 - 8 * y) * EXP(- y);
}
double D4Gauss(double x) {
  double y = x * x; 
  return ((16 * y - 48) * y + 12) * EXP(- y);
}




/* Whittle-Matern or Whittle or Besset ---- rescaled form of Whittle-Matern,
    see also there */ 

#define LOW_MATERN 1e-20
double logWM2(double x, double nu1, double nu2, double factor,
	      whittle_work_type *work) {
  
#if defined compatibility_to_R_h
  if (work == NULL) work = &(KEYT()->whittle);
#else
  assert(work != NULL);
#endif
  // check calling functions, like hyperbolic and gneiting if any changings !!

  //  printf("%10g %10g %10g %10g\n", x, nu1, nu2, factor);

   double v, y, loggamma,
    nu = 0.5 * (nu1 + nu2),
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
		   scale = 1.0;
  if (factor!=0.0) scale = factor * SQRT(nuThres);
  bool simple = nu1 == nu2 || nu > MATERN_NU_THRES;
  double bk[MATERN_NU_THRES + 1U];

  if (x > LOW_MATERN && nu < RF_INF) {
    if (x == RF_INF) return RF_NEGINF;
    if (simple) {
      if (nuThres != work->nuOld) {
	work->nuOld = nuThres;
	work->loggamma_old = LOGGAMMAFN(nuThres);
      } 
      loggamma = work->loggamma_old;      
    } else {
      if (nu1 != work->nu1old) {
	work->nu1old = nu1;
	work->loggamma1old = LOGGAMMAFN(nu1);
      }
      if (nu2 != work->nu2old) {
	work->nu2old = nu2;
	work->loggamma2old = LOGGAMMAFN(nu2);
      }
      loggamma = 0.5 * (work->loggamma1old + work->loggamma2old);
    }
    
    y = x  * scale;
    v = LOG2 + nuThres * LOG(0.5 * y) - loggamma + 
      LOG(bessel_k_ex(y, nuThres, 2.0, bk)) - y;
  } else v = 0.0;
    
  if (nu > MATERN_NU_THRES) { // factor!=0.0 && 
    double w, 
      g = MATERN_NU_THRES / nu;
    y = x * factor / 2;
    w = logGauss(y);

    //if (nu>100) printf("nu=%10g %10e %10e %10e\n", nu, v, g, w);

    v = v * g + (1.0 - g) * w;
    if (nu1 != nu2) { // consistenz zw. nu1, nu2 und nuThres wiederherstellen
      v += LOGGAMMAFN(nu)- 0.5 * (LOGGAMMAFN(nu1) + LOGGAMMAFN(nu2)); //!nuThres
    }
    
    // if (!R_FINITE(v)) ERR0("non-finite value in the whittle-matern model -- value of 'nu' is much too large");

    //if (nu>100) printf("v=%10g \n", v);
  }

  return v;
}


double WM(double x, double nu, double factor,
	   whittle_work_type *work) {
  // check calling functions, like hyperbolic and gneiting if any changings !!
  return EXP(logWM2(x, nu, nu, factor, work));
}

double DWM(double x, double nu, double factor,
	   whittle_work_type *work) { 
#if defined compatibility_to_R_h
   if (work == NULL) work = &(KEYT()->whittle);
#else
  assert(work != NULL);
#endif
  double   y, v,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
		   scale = 1.0;
  if (factor != 0.0) scale = factor * SQRT(nuThres);
  double bk[MATERN_NU_THRES + 1U];
  
  if (x > LOW_MATERN && nu < RF_INF) {
    if (x == RF_INF) return 0.0;
  if (nuThres!=work->nuOld) {
      work->nuOld = nuThres;
      work->loggamma_old = LOGGAMMAFN(nuThres);
    }
    y = x * scale;  
    v = - 2.0 * EXP(nuThres * LOG(0.5 * y) - work->loggamma_old + 
			     LOG(bessel_k_ex(y, nuThres - 1.0, 2.0, bk)) - y);
  } else {
    v = (nuThres > 0.5) ? 0.0 : (nuThres < 0.5) ? INFTY : 1.253314137;
  }
  v *= scale;

  if (nu > MATERN_NU_THRES) {
    double w, 
      g = MATERN_NU_THRES / nu;
    scale = factor / 2.0;
    y = x * scale;
    w = DGauss(y) * scale;
    v = v * g + (1.0 - g) * w;
  }
  return v;
}

double DDWM(double x, double nu, double factor,
	   whittle_work_type *work) { 
#if defined compatibility_to_R_h
  if (work == NULL) work = &(KEYT()->whittle);
#else
  assert(work != NULL);
#endif
  double  y, v,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
		   scale = 1.0;
  if (factor != 0.0) scale = factor * SQRT(nuThres);
  double scaleSq  = scale * scale,
		   bk[MATERN_NU_THRES + 1U];
  
  if (x > LOW_MATERN && nu < RF_INF) {
    if (x == RF_INF) return 0.0;
    if (nuThres!=work->nuOld) {
      work->nuAlt = nuThres;
      work->gamma = GAMMAFN(nuThres);
    }
    y = x * scale;
    v = POW(0.5 * y , nuThres - 1.0) / work->gamma *
      (- bessel_k_ex(y, nuThres - 1.0, 1.0, bk)
       + y * bessel_k_ex(y, nuThres - 2.0, 1.0, bk));
  } else {
    v = (nu > 1.0) ? -0.5 / (nu - 1.0) : INFTY;
  }
  v *= scaleSq;

  if (nu > MATERN_NU_THRES) {
    double w, 
      g = MATERN_NU_THRES / nu;
    scale = factor / 2.0;
    scaleSq = scale * scale;
    y = x * scale;
    w = DDGauss(y) * scaleSq;
    v = v * g + (1.0 - g) * w;
  }
  return v;
}

double D3WM(double x, double nu, double factor,
	   whittle_work_type *work) {
#if defined compatibility_to_R_h
  if (work == NULL) work = &(KEYT()->whittle);
#else
  assert(work != NULL);
#endif
  double y, v,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
    scale = (factor != 0.0) ? factor * SQRT(nuThres) : 1.0,
    scaleSq  = scale * scale;
   double bk[MATERN_NU_THRES + 1U];
 
  if (x > LOW_MATERN && nu < RF_INF) {
    if (x == RF_INF) return 0.0;
    if (nuThres!=work->nuOld) {
      work->nuAlt = nuThres;
      work->gamma = GAMMAFN(nuThres);
    }
    y = x * scale;
    v = POW(0.5 * y , nuThres - 1.0) / work->gamma *
      ( 3.0 * bessel_k_ex(y, nuThres - 2.0, 1.0, bk) 
	-y * bessel_k_ex(y, nuThres - 3.0, 1.0, bk)); 
  } else {
    v = 0.0;
  }
  v *= scaleSq * scale;
 
  if (nu > MATERN_NU_THRES) {
    double w, 
      g = MATERN_NU_THRES / nu;
    scale = factor / 2.0;
    scaleSq = scale * scale;
    y = x * scale;
    w = D3Gauss(y) * scaleSq * scale;
    v = v * g + (1.0 - g) * w;
  }
  return v;
}

double D4WM(double x,  double nu, double factor,
	   whittle_work_type *work) { 
#if defined compatibility_to_R_h
  if (work == NULL) work = &(KEYT()->whittle);
#else
  assert(work != NULL);
#endif
  double y, v,
    nuThres = nu < MATERN_NU_THRES ? nu : MATERN_NU_THRES,
    scale = (factor != 0.0) ? factor * SQRT(nuThres) : 1.0,
    scaleSq  = scale * scale;
  double bk[MATERN_NU_THRES + 1U];

  //  printf("x=%10g nu=%10g\n", x, nuThres);
  
  if (x > LOW_MATERN && nu < RF_INF) {
    if (x == RF_INF) return 0.0;
    if (nuThres!=work->nuOld) {
      work->nuAlt = nuThres;
      work->gamma = GAMMAFN(nuThres);
    }
    y = x * scale;
    v = 0.25 * POW(0.5 * y , nuThres - 3.0) / work->gamma *
      (+ 6.0 * (nuThres - 3.0 - y * y) * bessel_k_ex(y, nuThres - 3.0, 1.0, bk)
       + y * (3.0  + y * y) * bessel_k_ex(y, nuThres - 4.0, 1.0, bk)); 
  } else {
    v = INFTY;
    if (nuThres > 2.0) v = 0.75 / ((nuThres - 1.0) * (nuThres - 2.0));
  }
  v *= scaleSq * scaleSq;

  if (nu > MATERN_NU_THRES) {
    double w, 
      g = MATERN_NU_THRES / nu;
    scale = factor / 2.0;
    scaleSq = scale * scale;
    y = x * scale;
    w = D4Gauss(y) * scaleSq * scaleSq;
    v = v * g + (1.0 - g) * w;
  }

  // printf("v=%10g\n", v);
  
  return v;
}




/////////////////////////////////////
// not documente yet:




double incomplete_gamma(double start, double end, double s) {
  // int_start^end t^{s-1} e^{-t} \D t

  double
    v = 0.0, 
    w = 0.0;

  if (s <= 1.0) {
    if (start == 0.0) return RF_NAN;
  }
  
  double 
    e_start = EXP(-start),
    e_end = EXP(-end),
    power_start = POW(start, s),      
    power_end = end < RF_INF ? POW(end, s) : 0,
    factor = 1.0; 
  
  
  while (s < 0.0) {
    factor /= s;    
    v +=  factor * (power_end * e_end - power_start * e_start);
    power_start *= start;
    if (end < RF_INF) power_end *= end;
    s += 1.0;
  }

  w = pgamma(start, s, 1.0, 0, 0);  // q, shape, scale, lower, log


  if (R_FINITE(end)) w -= pgamma(end, s, 1.0, 0, 0);

  return v + GAMMAFN(s) * w * factor;
}

