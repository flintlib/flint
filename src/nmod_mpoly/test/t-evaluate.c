/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_evaluate, state)
{
    slong i, j, v;
    int tmul = 20;

    /* Check repeated evalone matches evalall */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g;
        mp_limb_t fe;
        mp_limb_t * vals;
        slong * perm;
        slong nvars, len;
        flint_bitcnt_t exp_bits;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);
        nvars = ctx->minfo->nvars;

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);

        len = n_randint(state, 50);
        exp_bits = n_randint(state, 200) + 1;

        perm = (slong *) flint_malloc(nvars*sizeof(slong));
        vals = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = n_randlimb(state);
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
            nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            fe = nmod_mpoly_evaluate_all_ui(f, vals, ctx);

            for (v = 0; v < nvars; v++)
            {
                nmod_mpoly_evaluate_one_ui(g, f, perm[v], vals[perm[v]], ctx);
                nmod_mpoly_assert_canonical(g, ctx);
                nmod_mpoly_evaluate_one_ui(f, f, perm[v], vals[perm[v]], ctx);
                nmod_mpoly_assert_canonical(f, ctx);
                if (!nmod_mpoly_equal(f, g, ctx))
                {
                    printf("FAIL\n");
                    flint_printf("Check evalone aliasing\ni: %wd  j: %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }
            if (!nmod_mpoly_equal_ui(f, fe, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check repeated evalone matches evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(vals);

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_ctx_clear(ctx);

        flint_free(perm);
    }

    /* Check add commutes with evalall */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, fg;
        mp_limb_t fe, ge, fge;
        mp_limb_t * vals;
        slong nvars, len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);
        nvars = ctx->minfo->nvars;

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(fg, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        vals = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = n_randlimb(state);
        }

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_add(fg, f, g, ctx);

            fe = nmod_mpoly_evaluate_all_ui(f, vals, ctx);
            ge = nmod_mpoly_evaluate_all_ui(g, vals, ctx);
            fge = nmod_mpoly_evaluate_all_ui(fg, vals, ctx);

            if (fge != nmod_add(fe, ge, ctx->mod))
            {
                printf("FAIL\n");
                flint_printf("Check add commutes with evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(vals);

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(fg, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check mul commutes with evalall */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, fg;
        mp_limb_t fe, ge, fge;
        mp_limb_t * vals;
        slong nvars, len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);
        nvars = ctx->minfo->nvars;

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(fg, ctx);

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        vals = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = n_randlimb(state);
        }

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_add(fg, f, g, ctx);

            fe = nmod_mpoly_evaluate_all_ui(f, vals, ctx);
            ge = nmod_mpoly_evaluate_all_ui(g, vals, ctx);
            fge = nmod_mpoly_evaluate_all_ui(fg, vals, ctx);

            if (fge != nmod_add(fe, ge, ctx->mod))
            {
                printf("FAIL\n");
                flint_printf("Check mul commutes with evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(vals);

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(fg, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
