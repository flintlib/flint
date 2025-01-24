/*
    Copyright (C) 2010,2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"

int
nmod_mat_solve_vec(nn_ptr x, const nmod_mat_t A, nn_srcptr b)
{
    nmod_mat_t X, B;
    int result;
    slong m;

    m = A->r;

    if (m == 0)
        return 1;

    B->entries = (nn_ptr) b;
    B->r = m;
    B->c = 1;
    B->stride = 1;
    B->mod = A->mod;

    X->entries = x;
    X->r = m;
    X->c = 1;
    X->stride = 1;
    X->mod = A->mod;

    result = nmod_mat_solve(X, A, B);

    return result;
}
