/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_siegel_randtest(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)
{
    slong g = arb_mat_nrows(tau);
    arb_mat_t re, im;
    slong k, j;

    arb_mat_init(re, g, g);
    arb_mat_init(im, g, g);

    for (k = 0; k < g; k++)
    {
        for (j = k; j < g; j++)
        {
            arb_randtest_precise(arb_mat_entry(re, k, j), state, prec, mag_bits);
            arb_set(arb_mat_entry(re, j, k), arb_mat_entry(re, k, j));
        }
    }

    arb_mat_randtest_spd(im, state, prec, mag_bits);
    acb_mat_set_real_imag(tau, re, im);

    arb_mat_clear(re);
    arb_mat_clear(im);
}
