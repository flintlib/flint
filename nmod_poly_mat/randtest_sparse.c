/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "nmod_poly.h"
#include "nmod_poly_mat.h"
#include "ulong_extras.h"

void
nmod_poly_mat_randtest_sparse(nmod_poly_mat_t A, flint_rand_t state, slong len,
    float density)
{
    slong i, j;

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            if (n_randint(state, 1000) < density * 1000)
            {
                slong l = n_randint(state, len + 1);
                l = FLINT_MAX(l, 1);
                nmod_poly_randtest(nmod_poly_mat_entry(A, i, j), state, l);
            }
            else
            {
                nmod_poly_zero(nmod_poly_mat_entry(A, i, j));
            }
        }
    }
}
