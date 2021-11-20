/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"


/* check total number of factors with multiplicity is between lower and upper */
void check_omega(slong lower, slong upper, const fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mod_mpoly_t q;
    fmpz_mod_mpoly_factor_t g, h;
    fmpz_t omega;

    fmpz_init(omega);
    fmpz_mod_mpoly_factor_init(g, ctx);
    fmpz_mod_mpoly_factor_init(h, ctx);
    fmpz_mod_mpoly_init(q, ctx);

    if (!fmpz_mod_mpoly_factor_wang(g, p, ctx))
    {
        flint_printf("FAIL:\nfactorization 1 could be computed\n");
        fflush(stdout);
        flint_abort();        
    }

    if (!fmpz_mod_mpoly_factor(h, p, ctx))
    {
        flint_printf("FAIL:\nfactorization 2 could be computed\n");
        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < g->num; i++)
    {
        if (g->poly[i].length < 1 || !fmpz_is_one(g->poly[i].coeffs + 0))
        {
            flint_printf("FAIL:\nfactorization is not unit normal\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, lower) < 0 || fmpz_cmp_si(omega, upper) > 0)
    {
        flint_printf("FAIL:\nfactorization has wrong number of factors\n");
        fflush(stdout);
        flint_abort();        
    }

    fmpz_mod_mpoly_factor_expand(q, g, ctx);
    if (!fmpz_mod_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL:\nfactorization does not match original polynomial\n");
        fflush(stdout);
        flint_abort();        
    }

    fmpz_mod_mpoly_factor_sort(g, ctx);
    fmpz_mod_mpoly_factor_sort(h, ctx);
    if (fmpz_mod_mpoly_factor_cmp(g, h, ctx) != 0)
    {
        flint_printf("FAIL:\nfactorizations do not match\n");
        fflush(stdout);
        flint_abort();        
    }

    for (i = 0; i < g->num; i++)
    {
        fmpz_mod_mpoly_factor(h, g->poly + i, ctx);
        if (h->num != 1 || !fmpz_is_one(h->exp + 0))
        {
            flint_printf("FAIL:\nfactor is reducible\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_mod_mpoly_clear(q, ctx);
    fmpz_mod_mpoly_factor_clear(g, ctx);
    fmpz_mod_mpoly_factor_clear(h, ctx);
    fmpz_clear(omega);
}


int
main(void)
{
    slong i, j, tmul = 30;
    FLINT_TEST_INIT(state);

    flint_printf("factor_wang....");
    fflush(stdout);

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong lower;
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, t;
        slong nfacs, len;
        ulong expbound, powbound, pow;
        slong n;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 7, 100);

        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        nfacs = 1 + (6 + n_randint(state, 6))/n;
        powbound = 1 + n_randint(state, 3);
        powbound = 1 + n_randint(state, powbound);
        expbound = 3 + 20/nfacs/n/powbound;

        lower = 0;
        fmpz_mod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            len = 1 + n_randint(state, 10/powbound);
            fmpz_mod_mpoly_randtest_bound(t, state, len, expbound, ctx);
            if (fmpz_mod_mpoly_is_zero(t, ctx))
                fmpz_mod_mpoly_one(t, ctx);
            pow = 1 + n_randint(state, powbound);
            if (!fmpz_mod_mpoly_is_fmpz(t, ctx))
                lower += pow;
            fmpz_mod_mpoly_pow_ui(t, t, pow, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
        }

        check_omega(lower, WORD_MAX, a, ctx);

        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
