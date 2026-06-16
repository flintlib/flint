/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "fmpz.h"
#include "gr.h"
#include "gr_generic.h"
#include "profiler.h"
#include "radix.h"
#include "radix_padic.h"
#include "padic.h"
#include "nmod.h"

int main()
{
    flint_rand_t state;
    ulong p;
    slong n;

    flint_rand_init(state);
    radix_t radix;

    p = 7;
    slong maxlen = 128;
    radix_init(radix, p, 0);

    flint_printf("precision %wu^n, limb radix = %wu^%wd\n\n", p, p, radix->exp);

    flint_printf("       n        len        naive       delayed\n\n");


    for (n = 10; n <= 11000000; n *= 2)
    {
        slong len;
        gr_ctx_t ctx;
        radix_padic_t r1, r2;
        radix_padic_struct *v1, *v2;
        fmpz_t fp, pn, fx;

        gr_ctx_init_radix_padic(ctx, p, n, RADIX_PADIC_PREC_INF, 0);

        fmpz_init_set_ui(fp, p);
        fmpz_init(pn);
        fmpz_init(fx);

        fmpz_ui_pow_ui(pn, p, n);

        v1 = gr_heap_init_vec(maxlen, ctx);
        v2 = gr_heap_init_vec(maxlen, ctx);
        radix_padic_init(r1, ctx);
        radix_padic_init(r2, ctx);

        for (slong i = 0; i < maxlen; i++)
        {
            fmpz_randm(fx, state, pn);
            if (0 && n_randint(state, 2))
                fmpz_neg(fx, fx);
            radix_padic_set_fmpz(v1 + i, fx, ctx);

            fmpz_randm(fx, state, pn);
            //fmpz_set_ui(fx, n_randint(state, 100));
            if (0 && n_randint(state, 2))
                fmpz_neg(fx, fx);
            radix_padic_set_fmpz(v2 + i, fx, ctx);

            //v2[i].v += n_randint(state, n / 4);
        }

        double t1, t2, FLINT_SET_BUT_UNUSED(__);

        for (len = 1; len <= maxlen; len *= 2)
        {
            TIMEIT_START;
            GR_MUST_SUCCEED(radix_padic_dot_strided_naive(r1, NULL, 0, v1, 1, v2, 1, len, ctx));
            TIMEIT_STOP_VALUES(__, t1);

            TIMEIT_START;
            GR_MUST_SUCCEED(radix_padic_dot_strided_delayed(r2, NULL, 0, v1, 1, v2, 1, len, ctx));
            TIMEIT_STOP_VALUES(__, t2);

            flint_printf("%8wd   %8wd   %10.2e    %10.2e   %5.2fx\n", n, len, t1, t2, t1 / t2);
        }

        radix_padic_clear(r1, ctx);
        radix_padic_clear(r2, ctx);
        gr_heap_clear_vec(v1, maxlen, ctx);
        gr_heap_clear_vec(v2, maxlen, ctx);

        fmpz_clear(fx);
        fmpz_clear(fp);
        fmpz_clear(pn);

        gr_ctx_clear(ctx);

        flint_printf("\n");
    }

    radix_clear(radix);
    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}

