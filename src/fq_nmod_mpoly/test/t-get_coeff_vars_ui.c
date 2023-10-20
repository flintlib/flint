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

TEST_FUNCTION_START(fq_nmod_mpoly_get_coeff_vars_ui, state)
{
    slong i, j1, j2;

    /* check simple example */
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g;
        /* get the coefficient of y^1*x^2*/
        slong varl[2] = {1, 0};
        ulong expl[2] = {1, 2};
        const char * vars[] = {"x", "y", "z", "w"};

        fq_nmod_mpoly_ctx_init_deg(ctx, 4, ORD_LEX, 11, 2);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_set_str_pretty(f, "x^2*y*(z+w)+x+y+x*y^2+z^2+w^2", vars, ctx);
        fq_nmod_mpoly_set_str_pretty(g, "z+w", vars, ctx);
        fq_nmod_mpoly_get_coeff_vars_ui(f, f, varl, expl, 2, ctx);
        if (!fq_nmod_mpoly_equal(f, g, ctx))
        {
            flint_printf("FAIL\ncheck simple example\n");
            fflush(stdout);
            flint_abort();
        }
        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* check 1 variable sum of coefficients */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, m;
        slong len, nvars;
        ulong exp_bound;
        slong vars[1];
        ulong exps[1];
        slong var1;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        nvars = ctx->minfo->nvars;
        if (nvars < 1)
        {
            fq_nmod_mpoly_ctx_clear(ctx);
            continue;
        }

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(m, ctx);

        len = n_randint(state, 200);
        exp_bound = n_randint(state, 20) + 1;

        fq_nmod_mpoly_randtest_bound(f, state, len, exp_bound, ctx);
        fq_nmod_mpoly_repack_bits(f, f, f->bits + n_randint(state, FLINT_BITS), ctx);

        var1 = n_randint(state, nvars);

        fq_nmod_mpoly_zero(h, ctx);
        for (j1 = 0; j1 < exp_bound; j1++)
        {
            vars[0] = var1;
            exps[0] = j1;
            fq_nmod_mpoly_get_coeff_vars_ui(g, f, vars, exps, 1, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);
            fq_nmod_mpoly_gen(m, var1, ctx);
            if (!fq_nmod_mpoly_pow_ui(m, m, j1, ctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check power success\ni = %wd, j1 = %wd\n", i, j1);
                fflush(stdout);
                flint_abort();
            }
            fq_nmod_mpoly_mul(g, g, m, ctx);
            fq_nmod_mpoly_add(h, h, g, ctx);
        }

        if (!fq_nmod_mpoly_equal(f, h, ctx))
        {
            flint_printf("FAIL\n");
            flint_printf("Check 1 variable sum of coefficients\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_clear(m, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* check 2 variable sum of coefficients */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, m;
        slong len, nvars;
        ulong exp_bound;
        slong vars[2];
        ulong exps[2];
        slong var1, var2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        nvars = ctx->minfo->nvars;
        if (nvars < 2)
        {
            fq_nmod_mpoly_ctx_clear(ctx);
            continue;
        }

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(m, ctx);

        len = n_randint(state, 200);
        exp_bound = n_randint(state, 12) + 1;

        fq_nmod_mpoly_randtest_bound(f, state, len, exp_bound, ctx);
        fq_nmod_mpoly_repack_bits(f, f, f->bits + n_randint(state, FLINT_BITS), ctx);

        var1 = n_randint(state, nvars - 1);
        var2 = 1 + var1 + n_randint(state, nvars - 1 - var1);

        fq_nmod_mpoly_zero(h, ctx);
        for (j1 = 0; j1 < exp_bound; j1++)
        {
        for (j2 = 0; j2 < exp_bound; j2++)
        {
            vars[0] = var1;
            exps[0] = j1;
            vars[1] = var2;
            exps[1] = j2;
            fq_nmod_mpoly_get_coeff_vars_ui(g, f, vars, exps, 2, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);
            fq_nmod_mpoly_gen(m, var1, ctx);
            if (!fq_nmod_mpoly_pow_ui(m, m, j1, ctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check power success\ni = %wd, j1 = %wd, j2 = %wd\n", i, j1, j2);
                fflush(stdout);
                flint_abort();
            }
            fq_nmod_mpoly_mul(g, g, m, ctx);
            fq_nmod_mpoly_gen(m, var2, ctx);
            if (!fq_nmod_mpoly_pow_ui(m, m, j2, ctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check power success\ni = %wd, j1 = %wd, j2 = %wd\n", i, j1, j2);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_mul(g, g, m, ctx);
            fq_nmod_mpoly_add(h, h, g, ctx);
        }
        }

        if (!fq_nmod_mpoly_equal(f, h, ctx))
        {
            flint_printf("FAIL\n"
                         "check 2 variable sum of coefficients\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_clear(m, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
