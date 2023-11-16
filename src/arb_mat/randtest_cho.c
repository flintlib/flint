/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"

void
arb_mat_randtest_cho(arb_mat_t mat, flint_rand_t state, slong prec, slong mag_bits)
{
    slong g = arb_mat_nrows(mat);
    slong i, j;

    arb_mat_zero(mat);
    for (i = 0; i < g; i++)
    {
        arb_randtest_positive(arb_mat_entry(mat, i, i), state, prec, mag_bits);
        for (j = 0; j < i; j++)
        {
            arb_randtest_precise(arb_mat_entry(mat, i, j), state, prec, mag_bits);
        }
    }
}
