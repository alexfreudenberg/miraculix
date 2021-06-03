
/*
 Authors 
 Alexander Freudenberg, afreuden@mail.uni-mannheim.de

Courtesy of Martin Schlather

 Copyright (C) 2021, Martin Schlather, Alexander Freudenberg

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, writne to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  
*/

#include <R_ext/Rdynload.h>
#include "rgpu.h"


#define CALLDEF_DO(name, n) {#name, (DL_FUNC) &name, n}
static R_CallMethodDef callMethods[] = {
   CALLDEF_DO(cu_matmul,3),
   CALLDEF_DO(access_pointer,1),
    {NULL, NULL, 0}
};

void R_init_rgpu(DllInfo *dll) {
    R_registerRoutines(dll,
            NULL,
            callMethods, // defined below
            NULL, // .Fortran
            NULL); // ext
    R_useDynamicSymbols(dll, FALSE);
}

void R_unload_rgpu(DllInfo *info) { }
