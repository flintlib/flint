/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly_factor.h"

void check_content(const fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, v;
    fq_nmod_mpoly_t q;
    fq_nmod_mpoly_factor_t g;
    fmpz_t deg;

    fq_nmod_mpoly_factor_init(g, ctx);
    fq_nmod_mpoly_init(q, ctx);
    fmpz_init(deg);

    if (!fq_nmod_mpoly_factor_content(g, p, ctx))
    {
        flint_printf("FAIL:\ncheck factorization could be computed\n");
        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < g->num; i++)
    {
        if (g->poly[i].length < 1 || g->poly[i].coeffs[0] != 1)
        {
            flint_printf("FAIL:\nfactorization is not unit normal\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fq_nmod_mpoly_factor_expand(q, g, ctx);
    if (!fq_nmod_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL:\nfactorization does not match original polynomial\n");
        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < g->num; i++)
    {
        for (v = 0; v < ctx->minfo->nvars; v++)
        {
            if (fq_nmod_mpoly_length(g->poly + i, ctx) < 2)
            {
                if (!fq_nmod_mpoly_is_gen(g->poly + i, -1, ctx))
                {
                    flint_printf("FAIL:\nmonomial is bad\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
            else
            {
                if (!fq_nmod_mpoly_content_vars(q, g->poly + i, &v, 1, ctx))
                {
                    flint_printf("FAIL:\ncheck content could be computed\n");
                    fflush(stdout);
                    flint_abort();
                }

                fq_nmod_mpoly_degree_fmpz(deg, g->poly + i, v, ctx);
                if (!fq_nmod_mpoly_is_one(q, ctx) && !fmpz_is_zero(deg))
                {
                    flint_printf("FAIL:\ncontent is bad\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }
    }

    fq_nmod_mpoly_clear(q, ctx);
    fq_nmod_mpoly_factor_clear(g, ctx);
    fmpz_clear(deg);
}

TEST_FUNCTION_START(fq_nmod_mpoly_factor_content, state)
{
    slong i, j, k, tmul = 30;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, t;
        slong n, nfacs, len;
        ulong * expbounds;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 5, FLINT_BITS, 4);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(t, ctx);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        nfacs = 5 + (8 + n_randint(state, 8))/n;
        expbounds = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, ulong);
        for (k = 0; k < ctx->minfo->nvars; k++)
            expbounds[k] = 3;

        fq_nmod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            len = 1 + n_randint(state, 9);
            if (ctx->minfo->nvars > 0)
            {
                k = n_randint(state, ctx->minfo->nvars);
                expbounds[k] = 1;
            }
            fq_nmod_mpoly_randtest_bounds(t, state, len, expbounds, ctx);
            if (!fq_nmod_mpoly_is_zero(t, ctx))
            {
                fq_nmod_mpoly_mul(t, a, t, ctx);
                if (t->length < 1500)
                    fq_nmod_mpoly_swap(a, t, ctx);
            }
        }

        check_content(a, ctx);

        flint_free(expbounds);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
