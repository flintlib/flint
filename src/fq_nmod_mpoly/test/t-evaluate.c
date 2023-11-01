/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_evaluate, state)
{
    slong i, j, v;
    slong tmul = 5;

    /* Check repeated evalone matches evalall */
    for (i = 0; i < 2*tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g;
        fq_nmod_t fe;
        fq_nmod_struct ** vals;
        slong * perm;
        slong nvars, len;
        flint_bitcnt_t exp_bits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        nvars = ctx->minfo->nvars;

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_init(fe, ctx->fqctx);

        len = n_randint(state, 50);
        exp_bits = n_randint(state, 200) + 1;

        perm = (slong *) flint_malloc(nvars*sizeof(slong));
        vals = (fq_nmod_struct **) flint_malloc(nvars*sizeof(fq_nmod_struct *));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = (fq_nmod_struct *) flint_malloc(sizeof(fq_nmod_struct));
            fq_nmod_init(vals[v], ctx->fqctx);
            fq_nmod_randtest(vals[v], state, ctx->fqctx);
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
            fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            fq_nmod_mpoly_evaluate_all_fq_nmod(fe, f, vals, ctx);

            for (v = 0; v < nvars; v++)
            {
                fq_nmod_mpoly_evaluate_one_fq_nmod(g, f, perm[v], vals[perm[v]], ctx);
                fq_nmod_mpoly_assert_canonical(g, ctx);
                fq_nmod_mpoly_evaluate_one_fq_nmod(f, f, perm[v], vals[perm[v]], ctx);
                fq_nmod_mpoly_assert_canonical(f, ctx);
                if (!fq_nmod_mpoly_equal(f, g, ctx))
                {
                    flint_printf("FAIL\n");
                    flint_printf("Check evalone aliasing\ni: %wd  j: %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }
            if (!fq_nmod_mpoly_equal_fq_nmod(f, fe, ctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check repeated evalone matches evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < ctx->minfo->nvars; v++)
        {
            fq_nmod_clear(vals[v], ctx->fqctx);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fq_nmod_clear(fe, ctx->fqctx);
        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);

        flint_free(perm);
    }

    /* Check add commutes with evalall */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, fg;
        fq_nmod_t fe, ge, te, fge;
        fq_nmod_struct ** vals;
        slong nvars, len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        nvars = ctx->minfo->nvars;

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(fg, ctx);
        fq_nmod_init(fe, ctx->fqctx);
        fq_nmod_init(ge, ctx->fqctx);
        fq_nmod_init(te, ctx->fqctx);
        fq_nmod_init(fge, ctx->fqctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        exp_bits1 = n_randint(state, 80) + 1;
        exp_bits2 = n_randint(state, 80) + 1;

        vals = (fq_nmod_struct **) flint_malloc(nvars*sizeof(fq_nmod_struct *));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = (fq_nmod_struct *) flint_malloc(sizeof(fq_nmod_struct));
            fq_nmod_init(vals[v], ctx->fqctx);
            fq_nmod_randtest(vals[v], state, ctx->fqctx);
        }

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_add(fg, f, g, ctx);

            fq_nmod_mpoly_evaluate_all_fq_nmod(fe, f, vals, ctx);
            fq_nmod_mpoly_evaluate_all_fq_nmod(ge, g, vals, ctx);
            fq_nmod_mpoly_evaluate_all_fq_nmod(fge, fg, vals, ctx);

            fq_nmod_add(te, fe, ge, ctx->fqctx);
            if (!fq_nmod_equal(fge, te, ctx->fqctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check add commutes with evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fq_nmod_clear(vals[v], ctx->fqctx);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fq_nmod_clear(fe, ctx->fqctx);
        fq_nmod_clear(ge, ctx->fqctx);
        fq_nmod_clear(te, ctx->fqctx);
        fq_nmod_clear(fge, ctx->fqctx);
        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(fg, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check mul commutes with evalall */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, fg;
        fq_nmod_t fe, ge, te, fge;
        fq_nmod_struct ** vals;
        slong nvars, len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        nvars = ctx->minfo->nvars;

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(fg, ctx);
        fq_nmod_init(fe, ctx->fqctx);
        fq_nmod_init(ge, ctx->fqctx);
        fq_nmod_init(te, ctx->fqctx);
        fq_nmod_init(fge, ctx->fqctx);

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);
        exp_bits1 = n_randint(state, 70) + 1;
        exp_bits2 = n_randint(state, 70) + 1;

        vals = (fq_nmod_struct **) flint_malloc(nvars*sizeof(fq_nmod_struct *));
        for (v = 0; v < nvars; v++)
        {
            vals[v] = (fq_nmod_struct *) flint_malloc(sizeof(fq_nmod_struct));
            fq_nmod_init(vals[v], ctx->fqctx);
            fq_nmod_randtest(vals[v], state, ctx->fqctx);
        }

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_mul(fg, f, g, ctx);

            fq_nmod_mpoly_evaluate_all_fq_nmod(fe, f, vals, ctx);
            fq_nmod_mpoly_evaluate_all_fq_nmod(ge, g, vals, ctx);
            fq_nmod_mpoly_evaluate_all_fq_nmod(fge, fg, vals, ctx);

            fq_nmod_mul(te, fe, ge, ctx->fqctx);
            if (!fq_nmod_equal(fge, te, ctx->fqctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check add commutes with evalall\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (v = 0; v < nvars; v++)
        {
            fq_nmod_clear(vals[v], ctx->fqctx);
            flint_free(vals[v]);
        }
        flint_free(vals);

        fq_nmod_clear(fe, ctx->fqctx);
        fq_nmod_clear(ge, ctx->fqctx);
        fq_nmod_clear(te, ctx->fqctx);
        fq_nmod_clear(fge, ctx->fqctx);
        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(fg, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
