/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

void
acb_mat_randtest(acb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits)
{
    slong i, j, density;

    density = n_randint(state, 100);

    if (n_randint(state, 2))
        for (i = 0; i < acb_mat_nrows(mat); i++)
            for (j = 0; j < acb_mat_ncols(mat); j++)
                if (n_randint(state, 100) < density)
                    acb_randtest(acb_mat_entry(mat, i, j), state, prec, mag_bits);
                else
                    acb_zero(acb_mat_entry(mat, i, j));
    else
        for (i = 0; i < acb_mat_nrows(mat); i++)
            for (j = 0; j < acb_mat_ncols(mat); j++)
                if (n_randint(state, 100) < density)
                    acb_randtest_precise(acb_mat_entry(mat, i, j), state, prec, mag_bits);
                else
                    acb_zero(acb_mat_entry(mat, i, j));
}

