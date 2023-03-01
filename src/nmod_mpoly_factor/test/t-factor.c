/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


/* check total number of factors with multiplicity is between lower and upper */
void check_omega(
    slong lower,
    slong upper,
    const nmod_mpoly_t p,
    const nmod_mpoly_ctx_t ctx,
    const char * s)
{
    slong i;
    nmod_mpoly_t q;
    nmod_mpoly_factor_t g, h;
    fmpz_t omega;

    fmpz_init(omega);
    nmod_mpoly_factor_init(g, ctx);
    nmod_mpoly_factor_init(h, ctx);
    nmod_mpoly_init(q, ctx);

    if (!nmod_mpoly_factor(g, p, ctx))
    {
        flint_printf("FAIL: %s\ncheck factorization could be computed\n", s);
        fflush(stdout);
        flint_abort();        
    }

    for (i = 0; i < g->num; i++)
    {
        if (g->poly[i].length < 1 || g->poly[i].coeffs[0] != 1)
        {
            flint_printf("FAIL: %s\nfactorization is not unit normal\n", s);
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, lower) < 0 || fmpz_cmp_si(omega, upper) > 0)
    {
        flint_printf("FAIL: %s\nfactorization has wrong number of factors\n", s);
        fflush(stdout);
        flint_abort();        
    }

    nmod_mpoly_factor_expand(q, g, ctx);
    if (!nmod_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL: %s\nfactorization does not match\n", s);
        fflush(stdout);
        flint_abort();        
    }

    for (i = 0; i < g->num; i++)
    {
        nmod_mpoly_factor(h, g->poly + i, ctx);
        if (h->num != 1 || !fmpz_is_one(h->exp + 0))
        {
            flint_printf("FAIL: %s\nfactor is reducible\n", s);
            fflush(stdout);
            flint_abort();
        }
    }

    nmod_mpoly_clear(q, ctx);
    nmod_mpoly_factor_clear(g, ctx);
    nmod_mpoly_factor_clear(h, ctx);
    fmpz_clear(omega);
}


int
main(void)
{
    slong i, j, tmul = 30;
    FLINT_TEST_INIT(state);

    flint_printf("factor....");
    fflush(stdout);

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a;
        const char * vars[] = {"x", "y", "z"};

        nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, 2);
        nmod_mpoly_init(a, ctx);

        nmod_mpoly_set_str_pretty(a,
                "x^4120+x^4118*y^2+x^3708*y^400+x^3706*y^402+x^2781*y^1300+"
                "x^2779*y^1302+x^1339*y^2700+x^927*y^3100+y^4000+x^7172*y^4167+"
                "x^8349*y^4432+x^8347*y^4434+x^6760*y^4567+x^5833*y^5467+"
                "x^5568*y^7132+x^11401*y^8599", vars, ctx);
        check_omega(2, 2, a, ctx, "irreducibility test");

        nmod_mpoly_set_str_pretty(a, "(1+x^2+y+x*y+y^2+z+x*z+y*z+z^2)*"
                "(1+x+x^2+x^3+y^2+x*y^2+y^3+y*z+x*y*z+z^2+x*z^2+y*z^2+z^3)*"
                "(1+x+x^4+x^5+y^3+x^2*y^3+y^5+y*z+x*y*z+x^2*y*z+x^3*y*z+y^3*z+"
                "x*y^3*z+y*z^2+x^2*y*z^2+z^3+x^2*z^3+y^2*z^3+z^5)", vars, ctx);
        check_omega(3, 3, a, ctx, "frobenius recombination");

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul*flint_test_multiplier(); i++)
    {
        slong lower;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, t;
        slong nfacs, len;
        ulong expbound, powbound, pow, expbounds[2];
        mp_limb_t p;

        p = n_randint(state, (i % 2 == 0) ? 4 : FLINT_BITS - 1) + 1;
        p = n_randbits(state, p);
        p = n_nextprime(p, 1);

        nmod_mpoly_ctx_init(ctx, 2, mpoly_ordering_randtest(state), p);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(t, ctx);

        nfacs = 5 + n_randint(state, 5);
        powbound = 1 + n_randint(state, 3);
        expbound = 3 + 70/nfacs/powbound;

        lower = 0;
        nmod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            len = 1 + n_randint(state, 10);
            expbounds[0] = 1 + n_randint(state, expbound);
            expbounds[1] = 1 + n_randint(state, expbound);
            nmod_mpoly_randtest_bounds(t, state, len, expbounds, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            pow = 1 + n_randint(state, powbound);
            if (!nmod_mpoly_is_ui(t, ctx))
                lower += pow;
            nmod_mpoly_pow_ui(t, t, pow, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
        }

        check_omega(lower, WORD_MAX, a, ctx, "bivariate");

        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong lower;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, t;
        slong n, nfacs, len;
        ulong expbound, powbound, pow;
        mp_limb_t p;

        p = n_randint(state, (i % 2 == 0) ? 4 : FLINT_BITS - 1) + 1;
        p = n_randbits(state, p);
        p = n_nextprime(p, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 10, p);

        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(t, ctx);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        nfacs = 1 + (6 + n_randint(state, 6))/n;
        powbound = 1 + n_randint(state, 3);
        powbound = 1 + n_randint(state, powbound);
        expbound = 3 + 100/nfacs/n/powbound;

        lower = 0;
        nmod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            len = 1 + n_randint(state, 10/powbound);
            nmod_mpoly_randtest_bound(t, state, len, expbound, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            pow = 1 + n_randint(state, powbound);
            if (!nmod_mpoly_is_ui(t, ctx))
                lower += pow;
            nmod_mpoly_pow_ui(t, t, pow, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
        }

        check_omega(lower, WORD_MAX, a, ctx, "multivariate");

        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
