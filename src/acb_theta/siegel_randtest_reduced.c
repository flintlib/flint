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
acb_siegel_randtest_reduced(acb_mat_t tau, flint_rand_t state, slong prec, slong mag_bits)
{
    slong g = acb_mat_nrows(tau);
    fmpz_mat_t mat;
    arb_t test;

    fmpz_mat_init(mat, 2 * g, 2 * g);
    arb_init(test);

    acb_siegel_randtest(tau, state, prec, mag_bits);
    acb_siegel_reduce(mat, tau, prec);
    acb_siegel_transform(tau, mat, tau, prec);

    fmpz_mat_clear(mat);
    arb_clear(test);
}
