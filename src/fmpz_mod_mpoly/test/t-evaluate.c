/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_evaluate, state)
{
    slong i, j, v;
    int tmul = 20;

    /* Check repeated evalone matches evalall */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g;
        fmpz_t fe;
        fmpz ** vals;
        slong * perm;
        slong nvars, len;
        flint_bitcnt_t exp_bits;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);
        nvars = ctx->minfo->nvars;

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_init(fe);

        len = n_randint(state, 50);
        exp_bits = n_randint(state, 200) + 1;

        perm = FLINT_ARRAY_ALLOC(nvars, slong);
        vals = FLINT_ARRAY_ALLOC(nvars, fmpz *);
        for (v = 0; v < nvars; v++)
        {
            vals[v] = FLINT_ARRAY_ALLOC(1, fmpz);
            perm[v] = v;
            fmpz_init(vals[v]);
            fmpz_randtest(vals[v], state, 200);
        }

        for (j = 0; j < 2*nvars; j++)
        {
            slong a = n_randint(state, nvars);
            slong b = n_randint(state, nvars);
            FLINT_SWAP(slong, perm[a], perm[b]);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_evaluate_all_fmpz(fe, f, vals, ctx);

            for (v = 0; v < nvars; v++)
            {
                fmpz_mod_mpoly_evaluate_one_fmpz(g, f, perm[v], vals[perm[v]], ctx);
                fmpz_mod_mpoly_assert_canonical(g, ctx);
                fmpz_mod_mpoly_evaluate_one_fmpz(f, f, perm[v], vals[perm[v]], ctx);
                fmpz_mod_mpoly_assert_canonical(f, ctx);
                if (!fmpz_mod_mpoly_equal(f, g, ctx))
                {
                    flint_printf("FAIL: Check evalone aliasing\n");
                    flint_printf("i: %wd  j: %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }

            if (!fmpz_mod_mpoly_equal_fmpz(f, fe, ctx))
            {
                flint_printf("FAIL: Check repeated evalone matches evalall\n");
                flint_printf("i: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fmpz_clear(vals[v]);
            flint_free(vals[v]);
        }
        flint_free(perm);
        flint_free(vals);

        fmpz_clear(fe);
        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);

    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, fg;
        fmpz_t fe, ge, fge;
        fmpz ** vals;
        slong nvars, len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);
        nvars = ctx->minfo->nvars;

        fmpz_init(fe);
        fmpz_init(ge);
        fmpz_init(fge);
        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(fg, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        exp_bits1 = n_randint(state, 100) + 1;
        exp_bits2 = n_randint(state, 100) + 1;

        vals = FLINT_ARRAY_ALLOC(nvars, fmpz *);
        for (v = 0; v < nvars; v++)
        {
            vals[v] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(vals[v]);
            fmpz_randtest(vals[v], state, 200);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fmpz_mod_mpoly_add(fg, f, g, ctx);

            fmpz_mod_mpoly_evaluate_all_fmpz(fe, f, vals, ctx);
            fmpz_mod_mpoly_evaluate_all_fmpz(ge, g, vals, ctx);
            fmpz_mod_mpoly_evaluate_all_fmpz(fge, fg, vals, ctx);

            fmpz_mod_add(fe, fe, ge, ctx->ffinfo);
            if (!fmpz_equal(fge, fe))
            {
                flint_printf("FAIL: Check add commutes with evalall\n");
                flint_printf("i: %wd  j: %wd\n", i, j);
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

        fmpz_clear(fe);
        fmpz_clear(ge);
        fmpz_clear(fge);
        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(fg, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, fg;
        fmpz_t fe, ge, fge;
        fmpz ** vals;
        slong nvars, len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);
        nvars = ctx->minfo->nvars;

        fmpz_init(fe);
        fmpz_init(ge);
        fmpz_init(fge);
        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(fg, ctx);

        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);
        exp_bits1 = n_randint(state, 100) + 1;
        exp_bits2 = n_randint(state, 100) + 1;

        vals = FLINT_ARRAY_ALLOC(nvars, fmpz *);
        for (v = 0; v < nvars; v++)
        {
            vals[v] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(vals[v]);
            fmpz_randtest(vals[v], state, 200);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fmpz_mod_mpoly_mul(fg, f, g, ctx);

            fmpz_mod_mpoly_evaluate_all_fmpz(fe, f, vals, ctx);
            fmpz_mod_mpoly_evaluate_all_fmpz(ge, g, vals, ctx);
            fmpz_mod_mpoly_evaluate_all_fmpz(fge, fg, vals, ctx);

            fmpz_mod_mul(fe, fe, ge, ctx->ffinfo);
            if (!fmpz_equal(fge, fe))
            {
                flint_printf("FAIL: Check mul commutes with evalall\n");
                flint_printf("i: %wd  j: %wd\n", i, j);
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

        fmpz_clear(fe);
        fmpz_clear(ge);
        fmpz_clear(fge);
        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(fg, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
