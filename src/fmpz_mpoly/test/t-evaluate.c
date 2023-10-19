/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_evaluate, state)
{
    slong i, j, v;

    {
        fmpz_t A1;
        fmpz_mpoly_t A, B;
        fmpz * Cp[3];
        fmpz C[3];
        fmpz_mpoly_ctx_t ctx;

        fmpz_mpoly_ctx_init(ctx, 3, ORD_LEX);

        fmpz_init(A1);
        fmpz_mpoly_init(B, ctx);
        fmpz_mpoly_init(A, ctx);
        for (i = 0; i < 3; i++)
        {
            Cp[i] = C + i;
            fmpz_init(C + i);
        }

        fmpz_mpoly_set_str_pretty(B,
                "1 + x1*x2^2 + x2^9999999999999999999999999*x3^9", NULL, ctx);

        fmpz_set_si(C + 0, 2);
        fmpz_set_si(C + 1, 2);
        fmpz_set_si(C + 2, 2);
        if (fmpz_mpoly_evaluate_all_fmpz(A1, B, Cp, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check non-example 1\n", i);
            fflush(stdout);
            flint_abort();
        }

        if (fmpz_mpoly_evaluate_one_fmpz(A, B, 1, C + 1, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check non-example 2\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_set_si(C + 0, 1);
        fmpz_set_si(C + 1, 1);
        fmpz_set_si(C + 2, 1);
        if (!fmpz_mpoly_evaluate_all_fmpz(A1, B, Cp, ctx) ||
            !fmpz_equal_si(A1, 3))
        {
            printf("FAIL\n");
            flint_printf("Check example 3\n", i);
            fflush(stdout);
            flint_abort();
        }

        if (!fmpz_mpoly_evaluate_one_fmpz(A, B, 1, C + 1, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check example 4\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_set_str_pretty(B, "1 + x1 + x3^9", NULL, ctx);
        if (!fmpz_mpoly_equal(A, B, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check example 4 equality\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(A1);
        fmpz_mpoly_clear(B, ctx);
        fmpz_mpoly_clear(A, ctx);
        for (i = 0; i < 3; i++)
            fmpz_clear(C + i);

        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check repeated evalone matches evalall */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;
        fmpz_t fe;
        fmpz ** vals;
        slong * perm;
        slong nvars, len1, exp_bound1;
        flint_bitcnt_t coeff_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = ctx->minfo->nvars;

        fmpz_mpoly_init(f, ctx);
        fmpz_init(fe);

        perm = (slong *) flint_malloc(nvars*sizeof(slong));

        len1 = n_randint(state, 50);
        exp_bound1 = n_randint(state, 10) + 1;
        coeff_bits = n_randint(state, 100) + 1;

        vals = (fmpz **) flint_malloc(nvars*sizeof(fmpz*));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(vals[v]);
            fmpz_randbits(vals[v], state, 10);
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
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);

            if (!fmpz_mpoly_evaluate_all_fmpz(fe, f, vals, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check evaluations success\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            for (v = 0; v < nvars; v++)
            {
                if (!fmpz_mpoly_evaluate_one_fmpz(f, f, perm[v], vals[perm[v]], ctx))
                {
                    printf("FAIL\n");
                    flint_printf("Check evaluations success\ni: %wd  j: %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
                fmpz_mpoly_assert_canonical(f, ctx);
            }
            if (!fmpz_mpoly_equal_fmpz(f, fe, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check repeated evalone matches evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fmpz_clear(vals[v]);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fmpz_mpoly_clear(f, ctx);

        fmpz_clear(fe);

        flint_free(perm);
    }

    /* Check multiprecision repeated evalone matches evalall */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;
        fmpz_t fe;
        fmpz ** vals;
        slong * perm;
        slong nvars, len1;
        flint_bitcnt_t exp_bits, coeff_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = ctx->minfo->nvars;

        fmpz_mpoly_init(f, ctx);
        fmpz_init(fe);

        perm = (slong *) flint_malloc(nvars*sizeof(slong));

        len1 = n_randint(state, 50);
        exp_bits = n_randint(state, 200) + 1;
        coeff_bits = n_randint(state, 200) + 1;

        vals = (fmpz **) flint_malloc(nvars*sizeof(fmpz*));
        for (v = 0; v < nvars; v++)
        {
            /* only evaluate at 0, 1, or -1 */
            vals[v] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(vals[v]);
            fmpz_set_si(vals[v], n_randint(state, UWORD(3)) - WORD(1));
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
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits, ctx);
            if (!fmpz_mpoly_evaluate_all_fmpz(fe, f, vals, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check evaluations success\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            for (v = 0; v < nvars; v++)
            {
                if (!fmpz_mpoly_evaluate_one_fmpz(f, f, perm[v], vals[perm[v]], ctx))
                {
                    printf("FAIL\n");
                    flint_printf("Check evaluations success\ni: %wd  j: %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }

                fmpz_mpoly_assert_canonical(f, ctx);
            }
            if (!fmpz_mpoly_equal_fmpz(f, fe, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check multiprecision repeated evalone matches evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fmpz_clear(vals[v]);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fmpz_mpoly_clear(f, ctx);

        fmpz_clear(fe);

        flint_free(perm);
    }

    /* Check addition commutes with evalall */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, fg;
        fmpz_t fe, ge, fge, t;
        fmpz ** vals;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits;
        slong n;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = ctx->minfo->nvars;

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(fg, ctx);

        fmpz_init(fe);
        fmpz_init(ge);
        fmpz_init(fge);
        fmpz_init(t);

        len1 = n_randint(state, 500);
        len2 = n_randint(state, 500);

        n = FLINT_MAX(WORD(1), nvars);
        exp_bound1 = n_randint(state, 5000/n/n) + 1;
        exp_bound2 = n_randint(state, 5000/n/n) + 1;

        coeff_bits = n_randint(state, 100);

        vals = (fmpz **) flint_malloc(nvars*sizeof(fmpz*));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(vals[v]);
            fmpz_randbits(vals[v], state, 10);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_add(fg, f, g, ctx);

            if (!fmpz_mpoly_evaluate_all_fmpz(fe, f, vals, ctx) ||
                !fmpz_mpoly_evaluate_all_fmpz(ge, g, vals, ctx) ||
                !fmpz_mpoly_evaluate_all_fmpz(fge, fg, vals, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check evaluations success\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_add(t, fe, ge);
            if (!fmpz_equal(t, fge))
            {
                printf("FAIL\n");
                flint_printf("Check addition commutes with evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fmpz_clear(vals[v]);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(fg, ctx);

        fmpz_clear(fe);
        fmpz_clear(ge);
        fmpz_clear(fge);
        fmpz_clear(t);

    }

    /* Check multiplication commutes with evalall */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, fg;
        fmpz_t fe, ge, fge, t;
        fmpz ** vals;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits;
        slong n;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = ctx->minfo->nvars;

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(fg, ctx);

        fmpz_init(fe);
        fmpz_init(ge);
        fmpz_init(fge);
        fmpz_init(t);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        n = FLINT_MAX(WORD(1), nvars);
        exp_bound1 = n_randint(state, 1000/n/n) + 1;
        exp_bound2 = n_randint(state, 1000/n/n) + 1;

        coeff_bits = n_randint(state, 100);

        vals = (fmpz **) flint_malloc(nvars*sizeof(fmpz*));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = (fmpz *) flint_malloc(sizeof(fmpz));
            fmpz_init(vals[v]);
            fmpz_randbits(vals[v], state, 10);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_mul_johnson(fg, f, g, ctx);

            if (!fmpz_mpoly_evaluate_all_fmpz(fe, f, vals, ctx) ||
                !fmpz_mpoly_evaluate_all_fmpz(ge, g, vals, ctx) ||
                !fmpz_mpoly_evaluate_all_fmpz(fge, fg, vals, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check evaluations success\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mul(t, fe, ge);
            if (!fmpz_equal(t, fge))
            {
                printf("FAIL\n");
                flint_printf("Check multiplication commutes with evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fmpz_clear(vals[v]);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(fg, ctx);

        fmpz_clear(fe);
        fmpz_clear(ge);
        fmpz_clear(fge);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
