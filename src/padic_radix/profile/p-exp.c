/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "double_extras.h"
#include "fmpz.h"
#include "gr.h"
#include "profiler.h"
#include "radix.h"
#include "padic_radix.h"
#include "padic.h"
#include "nmod.h"

int main()
{
    flint_rand_t state;
    ulong p;
    slong n;
    int pi;

    flint_rand_init(state);
    radix_t radix;

    for (pi = 0; pi < 2; pi++)
    {
        p = (pi == 0) ? 7 : n_nextprime(UWORD(1) << (FLINT_BITS - 1), 0);

        radix_init(radix, p, 0);

        flint_printf("precision %wu^n, limb radix = %wu^%wd\n\n", p, p, radix->exp);

        for (n = 10; n <= 11000000; n *= 2)
        {
            gr_ctx_t ctx;
            radix_integer_t x, y;
            padic_radix_t rx, ry;
            padic_ctx_t pctx;
            padic_t px, py;
            fmpz_t fp, pn, fx, fy;

            fmpz_init_set_ui(fp, p);
            fmpz_init(pn);
            fmpz_ui_pow_ui(pn, p, n);

            padic_ctx_init(pctx, fp, n - 1, n + 1, PADIC_SERIES);
            gr_ctx_init_padic_radix(ctx, p, n, PADIC_RADIX_PREC_INF, 0);

            radix_integer_init(x, radix);
            radix_integer_init(y, radix);

            padic_init2(px, n);
            padic_init2(py, n);

            fmpz_init(fx);
            fmpz_init(fy);

            padic_radix_init(rx, ctx);
            padic_radix_init(ry, ctx);

            do
            {
                fmpz_randm(fx, state, pn);
            } while (fmpz_divisible_ui(fx, p));

            padic_set_fmpz(px, fx, pctx);
            padic_shift(px, px, 1, pctx);

            padic_radix_set_fmpz(rx, fx, ctx);
            rx->N = n;
            rx->v = 1;
            rx->N += 1;

            double t1, t2 = 0.0, t3, t4, FLINT_SET_BUT_UNUSED(__);

            TIMEIT_START;
            padic_exp(py, px, pctx);
            TIMEIT_STOP_VALUES(__, t1);

            if (n < 1000)
            {
                TIMEIT_START;
                padic_radix_exp_rectangular(ry, rx, ctx);
                TIMEIT_STOP_VALUES(__, t2);
            }
            else
                t2 = D_NAN;

            TIMEIT_START;
            padic_radix_exp_balanced(ry, rx, ctx);
            TIMEIT_STOP_VALUES(__, t3);

            TIMEIT_START;
            padic_radix_exp(ry, rx, ctx);
            TIMEIT_STOP_VALUES(__, t4);

            if (n == 10)
                flint_printf("       n      padic     rectangular   balanced    default   speedup\n");

            flint_printf("%8wd   %10g  %10g  %10g  %10g  %6.3fx\n", n, t1, t2, t3, t4, t1 / t4);

            padic_radix_clear(rx, ctx);
            padic_radix_clear(ry, ctx);

            padic_clear(px);
            padic_clear(py);

            radix_integer_clear(x, radix);
            radix_integer_clear(y, radix);

            fmpz_clear(fx);
            fmpz_clear(fy);

            gr_ctx_clear(ctx);
            padic_ctx_clear(pctx);

            fmpz_clear(pn);
            fmpz_clear(fp);
        }

        radix_clear(radix);

        flint_printf("\n");
    }

    radix_clear(radix);
    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
