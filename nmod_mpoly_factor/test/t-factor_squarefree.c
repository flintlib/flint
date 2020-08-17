/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


void check_it(const nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx)
{
    slong i, j;
    nmod_mpoly_t q;
    nmod_mpoly_factor_t g, h;

    nmod_mpoly_factor_init(g, ctx);
    nmod_mpoly_factor_init(h, ctx);
    nmod_mpoly_init(q, ctx);

    if (!nmod_mpoly_factor_squarefree(g, p, ctx))
    {
        flint_printf("FAIL:\ncheck factorization could be computed\n");
        flint_abort();        
    }

    for (i = 0; i < g->num; i++)
    {
        if (g->poly[i].length < 1 || g->poly[i].coeffs[0] != 1)
        {
            flint_printf("FAIL:\nfactorization is not unit normal\n");
            flint_abort();
        }
    }

    nmod_mpoly_factor_expand(q, g, ctx);
    if (!nmod_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL:\nfactorization does not match original polynomial\n");
        flint_abort();        
    }

    for (i = 0; i < g->num; i++)
    {
        nmod_mpoly_factor_squarefree(h, g->poly + i, ctx);
        for (j = 0; j < h->num; j++)
        {
            if (!fmpz_is_one(h->exp + j))
            {
                flint_printf("FAIL:\nfactor has a square factor\n");
                flint_abort();
            }
        }
    }

    for (i = 1; i < g->num; i++)
    for (j = 0; j < i; j++)
    {
        if (!nmod_mpoly_gcd(q, g->poly + i, g->poly + j, ctx))
        {
            flint_printf("FAIL:\ncheck gcd could be computed\n");
        }

        if (!nmod_mpoly_is_one(q, ctx))
        {
            flint_printf("FAIL:\nbases have a common factor\n");
        }
    }

    nmod_mpoly_clear(q, ctx);
    nmod_mpoly_factor_clear(g, ctx);
    nmod_mpoly_factor_clear(h, ctx);
}


int
main(void)
{
    slong i, j, tmul = 30;
    FLINT_TEST_INIT(state);

    flint_printf("factor_squarefree....");
    fflush(stdout);

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, t;
        slong nfacs, len;
        ulong expbound, powbound, pow;
        mp_limb_t p;

        p = n_randint(state, (i % 2 == 0) ? 4 : FLINT_BITS - 1) + 1;
        p = n_randbits(state, p);
        p = n_nextprime(p, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 6, p);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(t, ctx);

        nfacs = 1 + (6 + n_randint(state, 6))/ctx->minfo->nvars;
        powbound = 1 + n_randint(state, 3);
        expbound = 3 + 25/nfacs/ctx->minfo->nvars/powbound;

        nmod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            len = 1 + n_randint(state, 60/powbound/ctx->minfo->nvars);
            nmod_mpoly_randtest_bound(t, state, len, expbound, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            pow = 1 + n_randint(state, powbound);
            nmod_mpoly_pow_ui(t, t, pow, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
        }

        check_it(a, ctx);

        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
