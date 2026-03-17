/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "mpoly.h"
#include "nmod_mpoly_factor.h"

/* Defined in t-factor.c, t-factor_wang.c, t-factor_zassenhaus.c and
 * t-factor_zippel.c */
#ifndef check_omega
#define check_omega check_omega
/* check total number of factors with multiplicity is between lower and upper */
void check_omega(
    slong lower,
    slong upper,
    const nmod_mpoly_t p,
    const nmod_mpoly_ctx_t ctx,
    int (* factor_fun)(nmod_mpoly_factor_t, const nmod_mpoly_t, const nmod_mpoly_ctx_t))
{
    slong i;
    nmod_mpoly_t q;
    nmod_mpoly_factor_t g, h;
    fmpz_t omega;

    fmpz_init(omega);
    nmod_mpoly_factor_init(g, ctx);
    nmod_mpoly_factor_init(h, ctx);
    nmod_mpoly_init(q, ctx);

    if (!factor_fun(g, p, ctx))
    {
        flint_printf("FAIL:\nfactorization 1 could be computed\n");
        fflush(stdout);
        flint_abort();
    }

    if (factor_fun != nmod_mpoly_factor)
    {
        if (!nmod_mpoly_factor(h, p, ctx))
        {
            flint_printf("FAIL:\nfactorization 2 could be computed\n");
            fflush(stdout);
            flint_abort();
        }
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

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, lower) < 0 || fmpz_cmp_si(omega, upper) > 0)
    {
        flint_printf("FAIL:\nfactorization has wrong number of factors\n");
        fflush(stdout);
        flint_abort();
    }

    nmod_mpoly_factor_expand(q, g, ctx);
    if (!nmod_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL:\nfactorization does not match original polynomial\n");
        fflush(stdout);
        flint_abort();
    }

    if (factor_fun != nmod_mpoly_factor)
    {
        nmod_mpoly_factor_sort(g, ctx);
        nmod_mpoly_factor_sort(h, ctx);
        if (nmod_mpoly_factor_cmp(g, h, ctx) != 0)
        {
            flint_printf("FAIL:\nfactorizations do not match\n");
            fflush(stdout);
            flint_abort();
        }
    }

    for (i = 0; i < g->num; i++)
    {
        nmod_mpoly_factor(h, g->poly + i, ctx);
        if (h->num != 1 || !fmpz_is_one(h->exp + 0))
        {
            flint_printf("FAIL:\nfactor is reducible\n");
            fflush(stdout);
            flint_abort();
        }
    }

    nmod_mpoly_clear(q, ctx);
    nmod_mpoly_factor_clear(g, ctx);
    nmod_mpoly_factor_clear(h, ctx);
    fmpz_clear(omega);
}
#endif

TEST_FUNCTION_START(nmod_mpoly_factor, state)
{
    slong i, j, tmul = 30;

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
        check_omega(2, 2, a, ctx, nmod_mpoly_factor); /* irreducibility test */

        nmod_mpoly_set_str_pretty(a, "(1+x^2+y+x*y+y^2+z+x*z+y*z+z^2)*"
                "(1+x+x^2+x^3+y^2+x*y^2+y^3+y*z+x*y*z+z^2+x*z^2+y*z^2+z^3)*"
                "(1+x+x^4+x^5+y^3+x^2*y^3+y^5+y*z+x*y*z+x^2*y*z+x^3*y*z+y^3*z+"
                "x*y^3*z+y*z^2+x^2*y*z^2+z^3+x^2*z^3+y^2*z^3+z^5)", vars, ctx);
        check_omega(3, 3, a, ctx, nmod_mpoly_factor); /* frobenius recombination */

        nmod_mpoly_set_str_pretty(a, "x^4*y^6*z^2+x^4*y^5+x^3*y*z^2+x*y^6*z^6+x*y^5*z^4+x*y^5*z^2+x*y^4+y*z^6+z^2", vars, ctx);
        check_omega(2, 2, a, ctx, nmod_mpoly_factor); /* deflate */

        nmod_mpoly_set_str_pretty(a, "(x^4*y^6*z^2+x^4*y^5+x^3*y*z^2+x*y^6*z^6+x*y^5*z^4+x*y^5*z^2+x*y^4+y*z^6+z^2)*(x^5*y^6*z^2+x^5*y^3*z^2+x^4*y^6*z^6+x^4*y^3*z^6+x^3*y^3*z^2+x^3*z^2+x^2*y^6+x*y^6*z^4+y^3)", vars, ctx);
        check_omega(4, 4, a, ctx, nmod_mpoly_factor); /* deflate */

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* issue 2611: crash in more: branch of n_bpoly_mod_factor_lgprime */
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a;
        const char * vars[] = {"x", "y", "z"};

        nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, 3);
        nmod_mpoly_init(a, ctx);

        nmod_mpoly_set_str_pretty(a,
                "x^18+x^12*y^6+x^6*y^12+y^18"
                "+x^12*y^4*z^2-x^10*y^6*z^2-x^6*y^10*z^2+x^4*y^12*z^2"
                "+x^12*y^2*z^4-x^10*y^4*z^4-x^4*y^10*z^4+x^2*y^12*z^4"
                "+x^12*z^6-x^10*y^2*z^6+x^6*y^6*z^6-x^2*y^10*z^6+y^12*z^6"
                "-x^6*y^2*z^10-x^4*y^4*z^10-x^2*y^6*z^10"
                "+x^6*z^12+x^4*y^2*z^12+x^2*y^4*z^12+y^6*z^12+z^18",
                vars, ctx);
        check_omega(1, 1, a, ctx, nmod_mpoly_factor);

        /* additional irreducible polys over GF(3) exercising the same path */
        nmod_mpoly_set_str_pretty(a,
                "x^11*y+2*x^10*y^2+2*x^10*y*z+2*x^10*z^2+x^9*y^3+x^9*y^2*z"
                "+x^9*z^3+2*x^8*y^4+x^8*y^3*z+x^8*y*z^3+2*x^8*z^4+x^7*y^5"
                "+x^7*y^4*z+x^7*y^3*z^2+2*x^7*y^2*z^3+x^7*y*z^4+2*x^7*z^5"
                "+x^6*y^6+x^6*y^5*z+2*x^6*y^4*z^2+2*x^6*y*z^5+x^5*y^6*z"
                "+2*x^5*y^5*z^2+2*x^5*y^4*z^3+2*x^5*y^2*z^5+2*x^5*z^7"
                "+x^4*y^8+x^4*y^7*z+2*x^4*y^6*z^2+x^4*y^5*z^3+2*x^4*y^3*z^5"
                "+x^4*y^2*z^6+2*x^4*z^8+x^3*y^9+x^3*y^8*z+2*x^3*y^6*z^3"
                "+2*x^3*y^5*z^4+x^3*y^4*z^5+x^3*y^3*z^6+2*x^3*y^2*z^7"
                "+2*x^3*y*z^8+2*x^3*z^9+2*x^2*y^7*z^3+x^2*y^6*z^4"
                "+2*x^2*y^5*z^5+x^2*y^4*z^6+x^2*y^3*z^7+x^2*y^2*z^8"
                "+x^2*y*z^9+x^2*z^10+x*y^10*z+x*y^9*z^2+x*y^8*z^3"
                "+x*y^6*z^5+x*y^4*z^7+x*y^2*z^9+2*y^12+2*y^8*z^4+y^7*z^5"
                "+y^6*z^6+y^5*z^7+2*y^2*z^10+y*z^11+2*z^12",
                vars, ctx);
        check_omega(1, 1, a, ctx, nmod_mpoly_factor);

        nmod_mpoly_set_str_pretty(a,
                "2*x^11*y+x^11*z+x^10*y*z+x^9*y^2*z+2*x^9*y*z^2+2*x^8*y^4"
                "+2*x^8*y^3*z+2*x^8*y^2*z^2+x^8*y*z^3+2*x^8*z^4+2*x^7*y^5"
                "+2*x^7*y^3*z^2+x^7*y^2*z^3+x^7*y*z^4+2*x^7*z^5+x^6*y^4*z^2"
                "+2*x^6*y^3*z^3+x^6*y^2*z^4+2*x^6*z^6+x^5*y^7+2*x^5*y^6*z"
                "+x^5*y^5*z^2+2*x^5*y^4*z^3+x^5*y^2*z^5+2*x^5*y*z^6+x^4*y^8"
                "+2*x^4*y^7*z+2*x^4*y^6*z^2+x^4*y^5*z^3+2*x^4*y^2*z^6"
                "+x^4*y*z^7+2*x^4*z^8+x^3*y^9+x^3*y^8*z+x^3*y^5*z^4"
                "+2*x^3*y^4*z^5+x^3*y^3*z^6+x^3*y*z^8+2*x^2*y^10+x^2*y^9*z"
                "+2*x^2*y^8*z^2+2*x^2*y^6*z^4+x^2*y^4*z^6+x^2*z^10+2*x*y^11"
                "+x*y^10*z+2*x*y^9*z^2+x*y^8*z^3+x*y^7*z^4+x*y^4*z^7"
                "+x*y^2*z^9+x*y*z^10+2*x*z^11+y^11*z+2*y^10*z^2+2*y^9*z^3"
                "+2*y^7*z^5+2*y^6*z^6+2*y^5*z^7+y^4*z^8+2*y^3*z^9+2*y*z^11",
                vars, ctx);
        check_omega(1, 1, a, ctx, nmod_mpoly_factor);

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
        ulong p;

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

        check_omega(lower, WORD_MAX, a, ctx, nmod_mpoly_factor); /* bivariate */

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
        ulong p;

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

        check_omega(lower, WORD_MAX, a, ctx, nmod_mpoly_factor); /* multivariate */

        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
