/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly_factor.h"

void check_it(const fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx)
{
    slong i, j;
    fmpz_mpoly_t q;
    fmpz_mpoly_factor_t g, h;

    fmpz_mpoly_factor_init(g, ctx);
    fmpz_mpoly_factor_init(h, ctx);
    fmpz_mpoly_init(q, ctx);

    if (!fmpz_mpoly_factor_squarefree(g, p, ctx))
    {
        flint_printf("FAIL:\ncheck factorization could be computed\n");
        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < g->num; i++)
    {
        if (g->poly[i].length < 1 || fmpz_sgn(g->poly[i].coeffs + 0) < 0)
        {
            flint_printf("FAIL:\nfactorization is not unit normal\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_mpoly_factor_expand(q, g, ctx);
    if (!fmpz_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL:\nfactorization does not match original polynomial\n");
        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < g->num; i++)
    {
        fmpz_mpoly_factor_squarefree(h, g->poly + i, ctx);
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
        if (!fmpz_mpoly_gcd(q, g->poly + i, g->poly + j, ctx))
        {
            flint_printf("FAIL:\ncheck gcd could be computed\n");
        }

        if (!fmpz_mpoly_is_one(q, ctx))
        {
            flint_printf("FAIL:\nbases have a common factor\n");
        }
    }

    fmpz_mpoly_clear(q, ctx);
    fmpz_mpoly_factor_clear(g, ctx);
    fmpz_mpoly_factor_clear(h, ctx);
}

TEST_FUNCTION_START(fmpz_mpoly_factor_squarefree, state)
{
    slong i, j, tmul = 30;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, t;
        flint_bitcnt_t coeff_bits;
        slong n, nfacs, len;
        ulong expbound, powbound, pow;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(t, ctx);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        nfacs = 1 + (7 + n_randint(state, 7))/n;
        powbound = 1 + n_randint(state, 5);
        expbound = 2 + 50/nfacs/n/powbound;

        fmpz_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            len = 1 + n_randint(state, 60/powbound/n);
            coeff_bits = 10 + n_randint(state, 200)/nfacs;
            fmpz_mpoly_randtest_bound(t, state, len, coeff_bits, expbound, ctx);
            if (fmpz_mpoly_is_zero(t, ctx))
                fmpz_mpoly_one(t, ctx);
            pow = 1 + n_randint(state, powbound);
            fmpz_mpoly_pow_ui(t, t, pow, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
        }

        check_it(a, ctx);

        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
