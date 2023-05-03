
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


#if !defined compatibility_to_R_h
#include <random>
#include <stdarg.h>
#endif


#include "def.h"
#include "Basic_RandomFieldsUtils.h"
#include "options_RFU.h"
#include "errors_messages.h"


#if !defined compatibility_to_R_h

double ownNA = 0.0,
  ownNaN = 0.0;


double gauss_random(double mu, double sigma) {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution<> d{mu,sigma};
  return (d(gen));
}

double uniform_random(double a, double b) {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::uniform_real_distribution<> d{a, b};
  return (d(gen));
}

double poisson_random(double lambda) {
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::poisson_distribution<> d{lambda};
  return (d(gen));
}

#endif


void stopIfNotIntI(Long i, Long line, const char *file) {
  if (i > MAXINT || i < -MAXINT)
    ERR3("value (%ld) not an integer at line %ld in %s\n", i, line, file);
}
  

void stopIfNotUIntI(Long i, Long line, const char *file) {
  if (i > UINT_MAX || i < 0)
    ERR3("value (%ld) not an unsigned integer at line %ld in %s\n", i, line, file);
}
  

void stopIfNotAnyIntI(Long i, Long line, const char *file) {
  if (i > INT_MAX || i < INT_MIN)
    ERR3("value (%ld) not an integer at line %ld in %s\n", i, line, file);
}
  
