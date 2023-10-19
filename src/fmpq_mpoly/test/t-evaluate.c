/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_evaluate, state)
{
    slong i, j, v;

    /* Check repeated evalone matches evalall */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        fmpq_t fe;
        fmpq ** vals;
        slong * perm;
        slong nvars;
        slong len1, exp_bound1;
        flint_bitcnt_t coeff_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = ctx->zctx->minfo->nvars;

        fmpq_mpoly_init(f, ctx);
        fmpq_init(fe);

        perm = (slong *) flint_malloc(nvars*sizeof(slong));

        len1 = n_randint(state, 50);
        exp_bound1 = n_randint(state, 10) + 1;
        coeff_bits = n_randint(state, 100) + 1;

        vals = (fmpq **) flint_malloc(nvars*sizeof(fmpq*));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = (fmpq *) flint_malloc(sizeof(fmpq));
            fmpq_init(vals[v]);
            fmpq_randtest(vals[v], state, 10);
            perm[v] = v;
        }

        for (j = 0; j < 2*nvars; j++)
        {
            slong a, b, c;
            a = n_randint(state, nvars);
            b = n_randint(state, nvars);
            c = perm[a];
            perm[a] = perm[b];
            perm[b] = c;
        }

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            if (!fmpq_mpoly_evaluate_all_fmpq(fe, f, vals, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check evaluation success\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            for (v = 0; v < nvars; v++)
            {
                if (!fmpq_mpoly_evaluate_one_fmpq(f, f, perm[v], vals[perm[v]], ctx))
                {
                    printf("FAIL\n");
                    flint_printf("Check evaluation success\ni: %wd  j: %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
                fmpq_mpoly_assert_canonical(f, ctx);
            }
            if (!fmpq_mpoly_equal_fmpq(f, fe, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check repeated evalone matches evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fmpq_clear(vals[v]);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fmpq_mpoly_clear(f, ctx);

        fmpq_clear(fe);

        flint_free(perm);
    }

    /* Check multiprecision repeated evalone matches evalall */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        fmpq_t fe;
        fmpq ** vals;
        slong * perm;
        slong nvars;
        slong len1;
        flint_bitcnt_t exp_bits, coeff_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = ctx->zctx->minfo->nvars;

        fmpq_mpoly_init(f, ctx);
        fmpq_init(fe);

        perm = (slong *) flint_malloc(nvars*sizeof(slong));

        len1 = n_randint(state, 50);
        exp_bits = n_randint(state, 200) + 1;
        coeff_bits = n_randint(state, 200) + 1;

        vals = (fmpq **) flint_malloc(nvars*sizeof(fmpq*));
        for (v = 0; v < nvars; v++)
        {
            /* only evaluate at 0, 1, or -1 */
            vals[v] = (fmpq *) flint_malloc(sizeof(fmpq));
            fmpq_init(vals[v]);
            fmpq_set_si(vals[v], n_randint(state, UWORD(3)) - WORD(1), UWORD(1));
            perm[v] = v;
        }

        for (j = 0; j < 2*nvars; j++)
        {
            slong a, b, c;
            a = n_randint(state, nvars);
            b = n_randint(state, nvars);
            c = perm[a];
            perm[a] = perm[b];
            perm[b] = c;
        }

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits, ctx);
            if (!fmpq_mpoly_evaluate_all_fmpq(fe, f, vals, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check evaluation success\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            for (v = 0; v < nvars; v++)
            {
                if (!fmpq_mpoly_evaluate_one_fmpq(f, f, perm[v], vals[perm[v]], ctx))
                {
                    printf("FAIL\n");
                    flint_printf("Check evaluation success\ni: %wd  j: %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }

                fmpq_mpoly_assert_canonical(f, ctx);
            }
            if (!fmpq_mpoly_equal_fmpq(f, fe, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check multiprecision repeated evalone matches evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fmpq_clear(vals[v]);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fmpq_mpoly_clear(f, ctx);

        fmpq_clear(fe);

        flint_free(perm);
    }

    /* Check addition commutes with evalall */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, fg;
        fmpq_t fe, ge, fge, t;
        fmpq ** vals;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        flint_bitcnt_t coeff_bits;
        slong n;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = ctx->zctx->minfo->nvars;

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(fg, ctx);

        fmpq_init(fe);
        fmpq_init(ge);
        fmpq_init(fge);
        fmpq_init(t);

        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        n = FLINT_MAX(WORD(1), nvars);
        exp_bound1 = n_randint(state, 2000/n/n) + 1;
        exp_bound2 = n_randint(state, 2000/n/n) + 1;

        coeff_bits = n_randint(state, 100);

        vals = (fmpq **) flint_malloc(nvars*sizeof(fmpq*));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = (fmpq *) flint_malloc(sizeof(fmpq));
            fmpq_init(vals[v]);
            fmpq_randbits(vals[v], state, 10);
        }

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpq_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpq_mpoly_add(fg, f, g, ctx);

            if (!fmpq_mpoly_evaluate_all_fmpq(fe, f, vals, ctx) ||
                !fmpq_mpoly_evaluate_all_fmpq(ge, g, vals, ctx) ||
                !fmpq_mpoly_evaluate_all_fmpq(fge, fg, vals, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check evaluation success\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpq_add(t, fe, ge);
            if (!fmpq_equal(t, fge))
            {
                printf("FAIL\n");
                flint_printf("Check addition commutes with evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fmpq_clear(vals[v]);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(fg, ctx);

        fmpq_clear(fe);
        fmpq_clear(ge);
        fmpq_clear(fge);
        fmpq_clear(t);
    }

    TEST_FUNCTION_END(state);
}
