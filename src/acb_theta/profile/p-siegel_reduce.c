/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "profiler.h"
#include "acb_mat.h"
#include "acb_theta.h"

static int usage(char * argv[])
{
    flint_printf("usage: %s g pstep pmax dstep dmax\n", argv[0]);
    return 1;
}

int main(int argc, char * argv[])
{
    slong g;
    slong prec, pmax, pstep;
    slong d, dmax, dstep;
    flint_rand_t state;
    acb_mat_t tau, w;
    arb_t r;
    fmpz_mat_t mat;
    slong j, k;

    if (argc < 6)
    {
        return usage(argv);
    }

    g = atol(argv[1]);
    pstep = atol(argv[2]);
    pmax = atol(argv[3]);
    dstep = atol(argv[4]);
    dmax = atol(argv[5]);

    flint_randinit(state);
    acb_mat_init(tau, g, g);
    acb_mat_init(w, g, g);
    arb_init(r);
    fmpz_mat_init(mat, 2 * g, 2 * g);

    acb_siegel_randtest_reduced(tau, state, pmax, 2);
    flint_printf("Starting matrix:\n");
    acb_mat_printd(tau, 5);

    for (prec = pstep; prec <= pmax; prec += pstep)
    {
        for (d = dstep; d <= dmax; d += dstep)
        {
            acb_mat_scalar_div_si(w, tau, d, prec);
            for (j = 0; j < g; j++)
            {
                for (k = 0; k <= j; k++)
                {
                    arb_urandom(r, state, prec);
                    acb_add_arb(acb_mat_entry(w, j, k), acb_mat_entry(w, j, k),
                        r, prec);
                    acb_set(acb_mat_entry(w, k, j), acb_mat_entry(w, j, k));
                }
            }

            flint_printf("prec = %wd, d = %wd\n", prec, d);

            TIMEIT_START;
            acb_siegel_reduce(mat, w, prec);
            TIMEIT_STOP;
        }
    }

    flint_randclear(state);
    acb_mat_clear(tau);
    acb_mat_clear(w);
    arb_clear(r);
    fmpz_mat_clear(mat);

    flint_cleanup();
    return 0;
}
