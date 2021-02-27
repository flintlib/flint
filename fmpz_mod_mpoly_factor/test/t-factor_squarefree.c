/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"


void check_it(const fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_mod_mpoly_t q;
    fmpz_mod_mpoly_factor_t g, h;

    fmpz_mod_mpoly_factor_init(g, ctx);
    fmpz_mod_mpoly_factor_init(h, ctx);
    fmpz_mod_mpoly_init(q, ctx);

    if (!fmpz_mod_mpoly_factor_squarefree(g, p, ctx))
    {
        flint_printf("FAIL:\ncheck factorization could be computed\n");
        flint_abort();        
    }

    for (i = 0; i < g->num; i++)
    {
        if (g->poly[i].length < 1 || !fmpz_is_one(g->poly[i].coeffs + 0))
        {
            flint_printf("FAIL:\nfactorization is not unit normal\n");
            flint_abort();
        }
    }

    fmpz_mod_mpoly_factor_expand(q, g, ctx);
    if (!fmpz_mod_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL:\nfactorization does not match original polynomial\n");
        flint_abort();
    }

    for (i = 0; i < g->num; i++)
    {
        fmpz_mod_mpoly_factor_squarefree(h, g->poly + i, ctx);
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
        if (!fmpz_mod_mpoly_gcd(q, g->poly + i, g->poly + j, ctx))
        {
            flint_printf("FAIL:\ncheck gcd could be computed\n");
        }

        if (!fmpz_mod_mpoly_is_one(q, ctx))
        {
            flint_printf("FAIL:\nbases have a common factor\n");
        }
    }

    fmpz_mod_mpoly_clear(q, ctx);
    fmpz_mod_mpoly_factor_clear(g, ctx);
    fmpz_mod_mpoly_factor_clear(h, ctx);
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
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, t;
        slong nfacs, len;
        ulong expbound, powbound, pow;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 6, 150);

        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        nfacs = 1 + (6 + n_randint(state, 6))/ctx->minfo->nvars;
        powbound = 1 + n_randint(state, 3);
        expbound = 3 + 25/nfacs/ctx->minfo->nvars/powbound;

        fmpz_mod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            len = 1 + n_randint(state, 80/powbound/ctx->minfo->nvars);
            fmpz_mod_mpoly_randtest_bound(t, state, len, expbound, ctx);
            if (fmpz_mod_mpoly_is_zero(t, ctx))
                fmpz_mod_mpoly_one(t, ctx);
            pow = 1 + n_randint(state, powbound);
            fmpz_mod_mpoly_pow_ui(t, t, pow, ctx);
            fmpz_mod_mpoly_mul(t, a, t, ctx);
            if (fmpz_mod_mpoly_length(t, ctx) < 500)
                fmpz_mod_mpoly_swap(a, t, ctx);
        }

        check_it(a, ctx);

        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
