
/*
 Authors 
 Alexander Freudenberg, afreuden@mail.uni-mannheim.de


 Copyright (C) 2021  Alexander Freudenberg

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


#ifndef MYPKG_H
#define MYPKG_H 1
#include <R.h>
#include <Rinternals.h>
#ifdef __cplusplus
extern "C" {
#endif
    SEXP cu_matmul(SEXP, SEXP, SEXP);
    SEXP access_pointer(SEXP);
#ifdef __cplusplus
}
#endif
#endif
