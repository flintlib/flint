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

/* check total number of factors with multiplicity is between lower and upper */
void check_it(const fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    fq_nmod_mpoly_t q;
    fq_nmod_mpoly_factor_t g, h;

    fq_nmod_mpoly_factor_init(g, ctx);
    fq_nmod_mpoly_factor_init(h, ctx);
    fq_nmod_mpoly_init(q, ctx);

    if (!fq_nmod_mpoly_factor_squarefree(g, p, ctx))
    {
        flint_printf("FAIL:\ncheck factorization 1 could be computed\n");
        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < g->num; i++)
    {
        if (!fq_nmod_mpoly_is_monic(g->poly + i, ctx))
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
        fq_nmod_mpoly_factor_squarefree(h, g->poly + i, ctx);
        for (j = 0; j < h->num; j++)
        {
            if (!fmpz_is_one(h->exp + j))
            {
                flint_printf("FAIL:\nfactor has a square factor\n");
                fflush(stdout);
                flint_abort();
            }
        }
    }

    for (i = 1; i < g->num; i++)
    for (j = 0; j < i; j++)
    {
        if (!fq_nmod_mpoly_gcd(q, g->poly + i, g->poly + j, ctx))
        {
            flint_printf("FAIL:\ncheck gcd could be computed\n");
        }

        if (!fq_nmod_mpoly_is_one(q, ctx))
        {
            flint_printf("FAIL:\nbases have a common factor\n");
        }
    }

    fq_nmod_mpoly_clear(q, ctx);
    fq_nmod_mpoly_factor_clear(g, ctx);
    fq_nmod_mpoly_factor_clear(h, ctx);
}

TEST_FUNCTION_START(fq_nmod_mpoly_factor_squarefree, state)
{
    slong i, j, tmul = 30;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, t;
        slong nfacs, len;
        ulong expbound, powbound, pow;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 5, FLINT_BITS, 4);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(t, ctx);

        nfacs = 2 + n_randint(state, 4);
        powbound = 3;
        expbound = 2 + 20/nfacs/FLINT_MAX(WORD(1), ctx->minfo->nvars);

        fq_nmod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            pow = 1 + n_randint(state, powbound);
            len = 1 + n_randint(state, 2 + 15/pow/nfacs);
            fq_nmod_mpoly_randtest_bound(t, state, len, expbound, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);

            fq_nmod_mpoly_pow_ui(t, t, pow, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
        }

        check_it(a, ctx);

        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
