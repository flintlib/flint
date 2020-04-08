/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by th e Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_sparse_vec.h"

void nmod_sparse_vec_window_init(nmod_sparse_vec_t window, const nmod_sparse_vec_t vec, slong i1, slong i2) 
{
    slong start, end;
    for (start = 0; start < vec->nnz && vec->entries[start].ind < i1; ++start);
    for (end = vec->nnz; end > 0 && vec->entries[end-1].ind >= i2; --end);
    window->entries = vec->entries + start;
    window->nnz = end - start;
}
