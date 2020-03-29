/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_mat.h"

void nmod_sparse_mat_set(nmod_sparse_mat_t mat, const nmod_sparse_mat_t src) 
{
    slong i;
    if(mat==src || mat->r==0) return;
    mat->r = src->r;
    mat->c = src->c;
    mat->c_off = src->c_off;
    mat->rows = realloc(mat->rows, mat->r*sizeof(*mat->rows));
    for(i=0; i<mat->r; ++i) nmod_sparse_vec_set(&mat->rows[i], &src->rows[i]);
}
