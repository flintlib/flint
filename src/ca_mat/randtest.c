/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

void
ca_mat_randtest(ca_mat_t mat, flint_rand_t state, slong length, slong bits, ca_ctx_t ctx)
{
    slong i, j, density;

    density = n_randint(state, 100);

    for (i = 0; i < ca_mat_nrows(mat); i++)
        for (j = 0; j < ca_mat_ncols(mat); j++)
            if (n_randint(state, 100) < density)
                ca_randtest(ca_mat_entry(mat, i, j), state, length, bits, ctx);
            else
                ca_zero(ca_mat_entry(mat, i, j), ctx);
}

void
ca_mat_randtest_rational(ca_mat_t mat, flint_rand_t state, slong bits, ca_ctx_t ctx)
{
    slong i, j, density;

    density = n_randint(state, 100);

    for (i = 0; i < ca_mat_nrows(mat); i++)
        for (j = 0; j < ca_mat_ncols(mat); j++)
            if (n_randint(state, 100) < density)
                ca_randtest_rational(ca_mat_entry(mat, i, j), state, bits, ctx);
            else
                ca_zero(ca_mat_entry(mat, i, j), ctx);
}
