/*
 Authors 
 Martin Schlather, schlather@math.uni-mannheim.de


 Copyright (C) 2018 -- 2018 Martin Schlather

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


#ifndef miraculix_utils_H
#define miraculix_utils_H 1

int *ToIntI(SEXP X, bool *create, bool round);
void freeGlobals();


#define ToInt(X)				\
  bool create##X = true;				\
  Uint *X##int = (Uint*) ToIntI(X, &create##X, false)

#define CondToInt(C,X)							\
  bool create##X = C;							\
  Uint *X##int = create##X ? (Uint*) ToIntI(X, &create##X, false) : NULL;

#endif
