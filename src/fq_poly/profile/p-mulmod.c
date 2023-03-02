/*
    Copyright 2010 William Hart
    Copyright 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "flint.h"
#include "ulong_extras.h"
#include "fq_poly.h"

typedef struct
{
    slong n;
    slong s;
    slong alg;
} info_t;

void
sample(void *arg, ulong count)
{
    info_t *info = (info_t *) arg;
    slong n = info->n, i, j, s = info->s, alg = info->alg;
    slong scale;

    fq_poly_t a, b, c, d, dinv;
    fq_ctx_t ctx;

    FLINT_TEST_INIT(state);
    
    fq_ctx_randtest(ctx, state);

    fq_poly_init2(a, n, ctx);
    fq_poly_init2(b, n, ctx);
    fq_poly_init2(c, 2 * n - 1, ctx);
    fq_poly_init2(d, s, ctx);
    fq_poly_init2(dinv, s, ctx);

    fq_poly_randtest_monic(a, state, n, ctx);
    fq_poly_randtest_monic(b, state, n, ctx);
    fq_poly_randtest_monic(d, state, s, ctx);

    fq_poly_reverse(dinv, d, s, ctx);
    fq_poly_inv_series_newton(dinv, dinv, s, ctx);

    scale = 1;
    if (n < 100000)
        scale = 10;
    if (n < 10000)
        scale = 100;
    if (n < 100)
        scale = 1000;

    for (i = 0; i < count; i++)
    {
        if (alg == 1)
        {
            prof_start();
            for (j = 0; j < scale; j++)
            {
                fq_poly_mulmod_preinv(c, a, b, d, dinv, ctx);
            }
            prof_stop();
        }
        else
        {
            prof_start();
            for (j = 0; j < scale; j++)
            {
                fq_poly_mulmod(c, a, b, d, ctx);
            }
            prof_stop();
        }
    }

    fq_poly_clear(a, ctx);
    fq_poly_clear(b, ctx);
    fq_poly_clear(c, ctx);
    fq_poly_clear(d, ctx);
    fq_poly_clear(dinv, ctx);
    fq_ctx_clear(ctx);
    FLINT_TEST_CLEANUP(state);
}

int
main(void)
{
    double min, max;
    info_t info;
    slong i, k, scale;


    for (k = 2; k <= 6; k++)
    {
        info.n = 1 << k;

        for (i = 0; i < 3; i++)
        {
            if (i == 0)
                info.s = info.n / 2;
            else if (i == 1)
                info.s = info.n + ((1 << (k + 1)) - (1 << k)) / 2 * i;
            else if (i == 2)
                info.s = info.n * 2 - 1;
            scale = 1;
            if (info.n < 100000)
                scale = 10;
            if (info.n < 10000)
                scale = 10;
            if (info.n < 100)
                scale = 10;

            info.alg = 1;
            prof_repeat(&min, &max, sample, (void *) &info);

            flint_printf
                ("length %wd, modulus degree %wd, min %.3g ms, max %.3g ms, norm %.3g\n",
                 info.n, info.s,
                 ((min / (double) FLINT_CLOCK_SCALE_FACTOR) / scale) /
                 2400000.0,
                 ((max / (double) FLINT_CLOCK_SCALE_FACTOR) / scale) /
                 2400000.0,
                 (((min / (double) FLINT_CLOCK_SCALE_FACTOR) / scale) /
                  2400000.0) * 500000.0 / info.n / FLINT_BIT_COUNT(info.n));
            fflush(stdout);

            info.alg = 2;
            prof_repeat(&min, &max, sample, (void *) &info);

            flint_printf
                ("length %wd, modulus degree %wd, min %.3g ms, max %.3g ms, norm %.3g\n",
                 info.n, info.s,
                 ((min / (double) FLINT_CLOCK_SCALE_FACTOR) / scale) /
                 2400000.0,
                 ((max / (double) FLINT_CLOCK_SCALE_FACTOR) / scale) /
                 2400000.0,
                 (((min / (double) FLINT_CLOCK_SCALE_FACTOR) / scale) /
                  2400000.0) * 500000.0 / info.n / FLINT_BIT_COUNT(info.n));
            fflush(stdout);
        }
    }

    return 0;
}
