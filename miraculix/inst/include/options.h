
/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2016 -- 2019 Martin Schlather

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


#ifndef miraculix_options_H
#define miraculix_options_H 1

#include <Rinternals.h>
#include "local1.h"
#include "local3.h"
#include "AutoMiraculix.h"

#define geneticsN 6
struct genetics_options {
  bool normalized, returnsigma, efficient;
  usr_bool centered;
  double digits;
  int method; // Alex: besserer Namen?
  // internal:
  double *pcentered;
  int ncentered;
};

// Alex: is 4 efficient?
#define genetics_START {			\
    true, false, true,					\
      True,						\
      3.0,						\
      (int) getAutoCodingIntern(),			\
      NULL, 0 /* internal */				\
      }

#define messagesN 1
extern const char * messages[messagesN];
struct messages_options{ 
  bool warn_address;
};
#define messages_START {			\
    true					\
      }

struct option_type {
  genetics_options genetics;
  messages_options messages;
};
extern option_type OPTIONS;

#define prefixN 2
extern const char * prefixlist[prefixN], **all[prefixN];
extern int allN[prefixN];

snpcoding getAutoCodingIntern();


typedef
struct KEY_type KEY_type;
struct KEY_type {
  KEY_type *next;
  option_type global;
  utilsoption_type global_utils;
  int pid,  visitingpid;
  bool ok, doshow;
  errorloc_type error_location;

  int *ToIntDummy;
  int  ToIntN = 0;
};
extern KEY_type *PIDKEY[PIDMODULUS];
KEY_type *KEYT();

  
#endif
  
  //option_type *global = &(KEYT()->global);
