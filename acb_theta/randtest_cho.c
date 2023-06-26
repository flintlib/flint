/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
arb_mat_randtest_cho(arb_mat_t mat, flint_rand_t state, slong prec,
                     slong mag_bits)
{
    slong g = arb_mat_nrows(mat);
    slong k, j;

    arb_mat_zero(mat);
    for (k = 0; k < g; k++)
    {
        arb_randtest_pos(arb_mat_entry(mat, k, k), state, prec, mag_bits);
        for (j = k + 1; j < g; j++)
        {
            arb_randtest_precise(arb_mat_entry(mat, k, j),
                                 state, prec, mag_bits);
        }
    }
}
