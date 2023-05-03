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


#ifdef TIME_AVAILABLE
#include <time.h>
#include <sys/time.h>
#endif


#ifdef TIME_AVAILABLE
#define MY_CLOCKS_PER_SEC CLOCKS_PER_SEC
// #define MY_CLOCKS_PER_SEC 750000
#define CLOCK_IN_SEC(X) ((double) (X) / (double) MY_CLOCKS_PER_SEC)
#define TIME_IN_SEC(X) ( (double) (X).tv_sec + (double) (X).tv_usec / 1000000.0 )
#define CLOCK(X) {							\
    clock_end = clock();						\
    gettimeofday(&TofD, NULL);	time_end = TIME_IN_SEC(TofD);		\
    PRINTF("%7.3f (%7.3f) sec for %s\n",				\
	   time_end-time_inter, CLOCK_IN_SEC(clock_end - clock_inter), X); \
    clock_inter = clock_end;						\
    time_inter = time_end;						\
  }
#define STARTCLOCK \
  clock_t clock_end, clock_start = clock(), clock_inter = clock_start;	\
  timeval TofD;								\
  gettimeofday(&TofD, NULL);						\
  double time_end, time_start = TIME_IN_SEC(TofD), time_inter = time_start; \

#define TOTALCLOCK(X) {							\
  CLOCK(X);								\
  PRINTF("%7.3f (%7.3f) sec total time\n",				\
	 time_end - time_start,  CLOCK_IN_SEC(clock_end - clock_start)); \
  }
#elif defined XXXXX99999
#define CLOCK(X) PRINTF("%s\n", X);
#define STARTCLOCK PRINT("startclock in %s, %d\n", __FILE__, __LINE__);
#define TOTALCLOCK PRINT("totalclock in %s, %d\n", __FILE__, __LINE__);
#else
#define CLOCK(X)
#define STARTCLOCK
#define TOTALCLOCK
#endif

#define CLOCK1(X,Y) { char msg[100]; SPRINTF(msg, X, Y); CLOCK(msg); }
#define CLOCK2(X,Y,Z) { char msg[100]; SPRINTF(msg, X, Y, Z); CLOCK(msg); }
#define CLOCK3(X,Y,Z,A) { char msg[100]; SPRINTF(msg, X, Y, Z, A); CLOCK(msg); }
 
