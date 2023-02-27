/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "mpoly.h"
#include "fmpz_mpoly.h"
#include "profiler.h"

int
main(void)
{
    slong i, j, k, total_time;
    timeit_t timer;
    flint_rand_t state;

    flint_randinit(state);

    flint_printf("------------------------------\n");
    total_time = 0;
    timeit_start(timer);
    for (i = 0; i < 20000; i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len1, len2;
        flint_bitcnt_t bits1, bits2;

        if (0 == i%1000)
        {
            flint_printf("%wd  ", i);
            fflush(stdout);
        }

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len1 = n_randint(state, 10) + 2;
        len2 = n_randint(state, 10) + 2;
        bits1 = n_randint(state, 9) + 1;
        bits2 = n_randint(state, 9) + 1;

        fmpz_mpoly_randtest_bits(h, state, len1, 10, bits1, ctx);
        fmpz_mpoly_randtest_bits(g, state, len2, 10, bits2, ctx);
        fmpz_mpoly_mul(f, h, g, ctx);

        if (mpoly_test_irreducible(f->exps, f->bits, f->length, ctx->minfo) &&
            h->length > 1 && g->length > 1)
        {
            flint_printf("FAIL: check reducible input\n");
            flint_printf("f: "); fmpz_mpoly_print_pretty(f, NULL, ctx); flint_printf("\n");
            flint_printf("g: "); fmpz_mpoly_print_pretty(g, NULL, ctx); flint_printf("\n");
            flint_printf("h: "); fmpz_mpoly_print_pretty(h, NULL, ctx); flint_printf("\n");
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }
    timeit_stop(timer);
    total_time += timer->wall;
    flint_printf("\nreducible time: %wd ms\n\n", timer->wall);

    flint_printf("------------------------------\n");
    total_time = 0;
    for (i = 2; i <= 16; i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;
        slong tot = 0, pos = 0;
        ulong * bounds = FLINT_ARRAY_ALLOC(i, ulong);

        fmpz_mpoly_ctx_init(ctx, i, ORD_LEX);
        fmpz_mpoly_init(f, ctx);

        timeit_start(timer);
        for (j = 0; j < 1000; j++)
        {
            flint_bitcnt_t bits = n_randint(state, 5) + 3;
            slong len = 10;

            for (k = 0; k < i; k++)
            {
                bounds[k] = n_urandint(state, UWORD(1) << bits) + 2;
                len += 3 + n_urandint(state, 2*bounds[k]);
            }

            fmpz_mpoly_randtest_bounds(f, state, len, 10, bounds, ctx);

            tot += 1;
            pos += mpoly_test_irreducible(f->exps, f->bits, f->length, ctx->minfo);
        }
        timeit_stop(timer);
        total_time += timer->wall;

        flint_printf("%wd vars: %f percent  %wd ms\n", i,
                                 (double)(pos)/(double)(tot)*100, timer->wall);

        flint_free(bounds);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }
    flint_printf("irreducible time: %wd ms\n\n", total_time);

    flint_randclear(state);

    return 0;
}

