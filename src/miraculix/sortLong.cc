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


#include "Basic_RFUlocal.h" // must be before anything else
#include "zzz_RFU.h"


typedef bool (*vergleich)(Long, Long, void *O);

bool smaller1L(Long i, Long j, void *orderd) {
  return ((double *) orderd)[i] < ((double *) orderd)[j];
}

bool greater1L(Long i, Long j, void *orderd) {
  return ((double *) orderd)[i] > ((double *) orderd)[j];
}

bool smallerLong1(Long i, Long j, void *orderedLong) {
  return ((Long *) orderedLong)[i] < ((Long *)orderedLong)[j];
}

bool greaterLong1(Long i, Long j, void *orderedLong) {
  return ((Long *)orderedLong)[i] > ((Long *)orderedLong)[j];
}


typedef bool (*vergleichX)(Long, Long, Long, void *O);
// vergleichX SMALLERXLong=NULL, GREATERXLong=NULL;

bool smallerL(Long i, Long j, Long orderDim, void *O)
{
  double *x, *y, *orderd = (double*) O;
  x = orderd + i * orderDim;
  y = orderd + j * orderDim;
  for(Long d=0; d<orderDim; d++)
     if (x[d] != y[d]) return x[d] < y[d];
  return false;
}

bool greaterL(Long i, Long j, Long orderDim, void *O)
{
  double *x, *y, *orderd = (double*) O;
  x = orderd + i * orderDim;
  y = orderd + j * orderDim;
  for(Long d=0; d<orderDim; d++)
    if (x[d] != y[d]) return x[d] > y[d];
  return false;
}


bool smallerLong(Long i, Long j, Long orderDim, void *O)
{
  Long *x, *y, *orderedLong = (Long*) O;
  x = orderedLong + i * orderDim;
  y = orderedLong + j * orderDim;
  for(Long d=0; d<orderDim; d++) {
    if (x[d] != y[d]) {
     return x[d] < y[d];
    }
  }
  return false;
}

bool greaterLong(Long i, Long j, Long orderDim, void *O)
{
  Long *x, *y, *orderedLong = (Long*) O;
  x = orderedLong + i * orderDim;
  y = orderedLong + j * orderDim;
  for(Long d=0; d<orderDim; d++)
    if (x[d] != y[d]) return x[d] > y[d];
  return false;
}


void orderLong(Long *pos, Long start, Long end,
	    vergleich SMALLER, vergleich GREATER,
	    void * orderd, Long order_from, Long order_to) {
  Long randpos, pivot, left, right, pivotpos;
  
  if( start < end ) {   
    randpos = (start + end) / 2;
    pivot = pos[randpos];
    pos[randpos] = pos[start];
    pos[start] = pivot;
    
    pivotpos=start; 
    left = start;
    right=end+1;   
    while (left < right) {      
      while (++left < right && SMALLER(pos[left], pivot, orderd)) pivotpos++;
      while (--right > left && GREATER(pos[right], pivot, orderd));      
      if (left < right) {
	Long swap=pos[left]; pos[left]=pos[right]; pos[right]=swap;
	pivotpos++;
      }
    }
    pos[start] = pos[pivotpos];
    pos[pivotpos] = pivot;
    if (start <= order_to && pivotpos > order_from)
      orderLong(pos, start, pivotpos-1, SMALLER, GREATER,
	    orderd,  order_from,  order_to);
    if (pivotpos < order_to && end >= order_from)
      orderLong(pos, pivotpos + 1, end, SMALLER, GREATER,
	    orderd, order_from, order_to);
  }
}


void XorderLong(Long *pos, Long start, Long end,
	    vergleichX SMALLER, vergleichX GREATER,
	    Long D, void * orderd, Long order_from, Long order_to ) {
  Long randpos, pivot, left, right, pivotpos;
  
  if( start < end ) {   
    randpos = (start + end) / 2;
    pivot = pos[randpos];
    pos[randpos] = pos[start];
    pos[start] = pivot;
    
    pivotpos=start; 
    left = start;
    right=end+1;   
    while (left < right) {
      while (++left < right && SMALLER(pos[left], pivot, D, orderd)) pivotpos++;
      while (--right > left && GREATER(pos[right], pivot, D, orderd));
      if (left < right) {
	Long swap=pos[left]; pos[left]=pos[right]; pos[right]=swap;
	pivotpos++;
      }
    }
    pos[start] = pos[pivotpos];
    pos[pivotpos] = pivot;
    if (start <= order_to && pivotpos > order_from)
      XorderLong(pos, start, pivotpos-1, SMALLER, GREATER,
	     D, orderd, order_from, order_to);
    if (pivotpos < order_to && end >= order_from)
      XorderLong(pos, pivotpos + 1, end, SMALLER, GREATER,
	     D, orderd, order_from, order_to);
  }
}


void orderingFromToL(double *d, Long len, int dim, Long *pos, Long from,
		     Long to,  usr_bool NAlast) 
{
  Long start, end;
  if (NAlast == Nan) {
    for (Long i=0; i<len; i++) pos[i]=i;
    end = len-1;
    start = 0;
  } else {
    if (dim != 1) ERR0("NAs only allowed for scalars");
    if (NAlast == True) {
       start = 0;
       end = -1;
       Long NAend = len;
       for (Long i=0; i<len; i++) 
	 if (ISNA(d[i]) || ISNAN(d[i])) pos[--NAend] = i;
	 else pos[++end] = i;
       assert(NAend - 1 == end);
    } else { // if (NAlast == False) {
      start = len;
      end = len -1;
      Long NAstart = -1;
      for (Long i=0; i<len; i++) 
	if (ISNA(d[i]) || ISNAN(d[i])) pos[++NAstart] = i;
	else pos[--start] = i;
      assert(NAstart + 1 == start);
    }
  }
  if (dim == 1) {
    orderLong(pos, start, end, smaller1L, greater1L,
	      (void *) d, from - 1, to - 1);
  } else {
    XorderLong(pos, start, end, smallerL, greaterL, dim,
	       (void*) d, from - 1, to - 1);
  }
}


void orderingL(double *d, Long len, int dim, Long *pos) {
  orderingFromToL(d, len, dim, pos, 1, len, Nan);
}


void orderingLongFromTo(Long *d, Long len, int dim, Long *pos, Long from,
			   Long to, usr_bool NAlast) {
  /* quicksort algorithm, slightly modified:
     does not sort the data, but d[pos] will be ordered 
     NOTE: pos must have the values 0,1,2,...,start-end !
     (orderdouble is a kind of sorting pos according to
     the variable d)
  */  
    Long start, end;
  if (NAlast == Nan) {
    for (Long i=0; i<len; i++) pos[i]=i;
    end = len-1;
    start = 0;
  } else {
    if (dim != 1) ERR0("NAs only allowed for scalars");
    if (NAlast == True) {
       start = 0;
       end = -1;
       Long NAend = len;
       for (Long i=0; i<len; i++) 
	 if (d[i] == NA_LONG) pos[--NAend] = i;
	 else pos[++end]=i;
       if (NAend - 1 != end) BUG;
    } else { // if (NAlast == False) {
      start = len;
      end = len -1;
      Long NAstart = -1;
      for (Long i=0; i<len; i++) 
	if (d[i] == NA_LONG) pos[++NAstart] = i;
	else pos[--start]=i;
      if (NAstart + 1 != start) BUG;
    }
  }
  if (dim == 1) {
    orderLong(pos, start, end, smallerLong1, greaterLong1, (void *) d,
	      from-1, to-1);
  } else {
    XorderLong(pos, start, end, smallerLong, greaterLong, dim, (void*) d,
	       from-1, to-1);
  }  
}
 
void orderingLong(Long *d, Long len, int dim, Long *pos)  {
  orderingLongFromTo(d, len, dim, pos, 1, len, Nan);
}




void quicksortL(Long start, Long end, double *orderd, Long order_from,
		   Long order_to)
{

  Long left, right, pivotpos;
  
  if( start < end ) {   
    Long randpos = (start + end) / 2;
    double pivot = orderd[randpos];
    orderd[randpos] = orderd[start];
    orderd[start] = pivot;
      
    pivotpos=start; 
    left = start;
    right = end+1;   
    
    while (left < right) {      
      while (++left < right && orderd[left] < pivot) pivotpos++;
      while (--right > left && orderd[right] > pivot);      
      if (left < right) {
	double swap = orderd[left];
	orderd[left]=orderd[right]; 
	orderd[right]=swap;
	pivotpos++;
      }
    }
    orderd[start] = orderd[pivotpos];
    orderd[pivotpos] = pivot;
    if (start <= order_to && pivotpos > order_from)
      quicksortL(start, pivotpos-1, orderd, order_from, order_to);
    if (pivotpos < order_to && end >= order_from)
      quicksortL(pivotpos + 1, end, orderd, order_from, order_to);
  }
}



void sortingFromToL(double *d, Long len, Long from, Long to, usr_bool NAlast) {
   Long start, end;
  if (NAlast == Nan) {
    end = len-1;
    start = 0;
  } if (NAlast == True) {
    start = end = 0;
    Long NAend = len - 1;
    while (end < NAend) {
      while (NAend >= 0 && (ISNA(d[NAend]) || ISNAN(d[NAend]))) NAend--;
      while (end < NAend && !ISNA(d[end]) && !ISNAN(d[end])) end++;
      if (end < NAend) {
	double swap = d[end];
	d[end] = d[NAend]; 
	d[NAend--] = swap;
      }
    }
     assert(NAend == end && false);
  } else { // if (NAlast == False) {
    start = end = len - 1;
    Long NAstart = 0;
    while (start > NAstart) {
      while(NAstart < len && (ISNA(d[NAstart]) || ISNAN(d[NAstart]))) NAstart++;
      while (start > NAstart && !ISNA(d[start]) && !ISNAN(d[start])) start--;
      if (start > NAstart) {
	double swap = d[start];
	d[start] = d[NAstart]; 
	d[NAstart++] = swap;
      }
    }   
   assert(NAstart == start);
  }
  quicksortL(start, end, d, from - 1, to - 1);
}

void sortingL(double *d, Long len, usr_bool NAlast) {
  sortingFromToL(d, len, 1, len, NAlast);
}

void sortLong(Long start, Long end, Long *orderedLong, Long order_from,
	      Long order_to) {

  Long left, right, pivotpos;
  
  if( start < end ) {   
    Long randpos = (start + end) / 2;
    Long pivot = orderedLong[randpos];
    orderedLong[randpos] = orderedLong[start];
    orderedLong[start] = pivot;
      
    pivotpos=start; 
    left = start;
    right = end+1;   
    
    while (left < right) {      
      while (++left < right && orderedLong[left] < pivot) pivotpos++;
      while (--right > left && orderedLong[right] > pivot);      
      if (left < right) {
	Long swap = orderedLong[left];
	orderedLong[left]=orderedLong[right]; 
	orderedLong[right]=swap;
	pivotpos++;
      }
    }
    orderedLong[start] = orderedLong[pivotpos];
    orderedLong[pivotpos] = pivot;
    if (start <= order_to && pivotpos > order_from)
      sortLong(start, pivotpos-1, orderedLong, order_from, order_to);
    if (pivotpos < order_to && end >= order_from)
      sortLong(pivotpos + 1, end, orderedLong, order_from, order_to);
  }
}



void sortingLongFromTo(Long *d, Long len, Long from, Long to,  usr_bool NAlast){
  /* quicksort algorithm, slightly modified:
     does not sort the data, but d[pos] will be ordered 
     NOTE: pos must have the values 0,1,2,...,start-end !
     (orderdouble is a kind of sorting pos according to
     the variable d)
  */  

   Long start, end;
  if (NAlast == Nan) {
    end = len-1;
    start = 0;
  } if (NAlast == True) {
    start = end = 0;
    Long NAend = len - 1;
    while (end < NAend) {
      while (NAend >= 0 && d[NAend] == NA_LONG) NAend--;
      while (end < NAend && d[end] != NA_LONG) end++;
      if (end < NAend) {
	Long swap = d[end];
	d[end] = d[NAend]; 
	d[NAend--] = swap;
      }
    }
     assert(NAend == end && false);
  } else { // if (NAlast == False) {
    start = end = len - 1;
    Long NAstart = 0;
    while (start > NAstart) {
      while(NAstart < len && d[NAstart] == NA_LONG) NAstart++;
      while (start > NAstart && d[start] != NA_LONG) start--;
      if (start > NAstart) {
	Long swap = d[start];
	d[start] = d[NAstart]; 
	d[NAstart++] = swap;
      }
    }  
    assert(NAstart == start);
  }
  sortLong(start, end, d, from - 1, to - 1);
}
 
void sortingLong(Long *d, Long len, usr_bool NAlast) {
  sortingLongFromTo(d, len, 1, len, NAlast);
}

