/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fq_nmod_mpoly_factor.h"

/* Defined in t-factor.c, t-factor_wang.c, t-factor_zassenhaus.c,
 * t-factor_zippel.c */
#ifndef check_omega
#define check_omega check_omega
/* check total number of factors with multiplicity is between lower and upper */
void check_omega(
        slong lower,
        slong upper,
        const fq_nmod_mpoly_t p,
        const fq_nmod_mpoly_ctx_t ctx,
        int (* factor_fun)(fq_nmod_mpoly_factor_t, const fq_nmod_mpoly_t, const fq_nmod_mpoly_ctx_t))
{
    slong i;
    fq_nmod_mpoly_t q;
    fq_nmod_mpoly_factor_t g, h;
    fmpz_t omega;

    fmpz_init(omega);
    fq_nmod_mpoly_factor_init(g, ctx);
    fq_nmod_mpoly_factor_init(h, ctx);
    fq_nmod_mpoly_init(q, ctx);

    if (!factor_fun(g, p, ctx))
    {
        flint_printf("FAIL:\ncheck factorization 1 could be computed\n");
        fflush(stdout);
        flint_abort();
    }

    if (factor_fun != fq_nmod_mpoly_factor)
    {
        if (!fq_nmod_mpoly_factor(h, p, ctx))
        {
            flint_printf("FAIL:\ncheck factorization 2 could be computed\n");
            fflush(stdout);
            flint_abort();
        }
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

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, lower) < 0 || fmpz_cmp_si(omega, upper) > 0)
    {
        flint_printf("FAIL:\nfactorization has wrong number of factors\n");
        fflush(stdout);
        flint_abort();
    }

    fq_nmod_mpoly_factor_expand(q, g, ctx);

    if (!fq_nmod_mpoly_equal(q, p, ctx))
    {
        flint_printf("FAIL:\nfactorization does not match original polynomial\n");
        fflush(stdout);
        flint_abort();
    }

    if (factor_fun != fq_nmod_mpoly_factor)
    {
        fq_nmod_mpoly_factor_sort(g, ctx);
        fq_nmod_mpoly_factor_sort(h, ctx);

        if (fq_nmod_mpoly_factor_cmp(g, h, ctx) != 0)
        {
            flint_printf("factorizations do not match\n");
            fflush(stdout);
            flint_abort();
        }
    }

    for (i = 0; i < g->num; i++)
    {
        fq_nmod_mpoly_factor(h, g->poly + i, ctx);
        if (h->num != 1 || !fmpz_is_one(h->exp + 0))
        {
            flint_printf("FAIL:\nfactor is reducible\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fq_nmod_mpoly_clear(q, ctx);
    fq_nmod_mpoly_factor_clear(g, ctx);
    fq_nmod_mpoly_factor_clear(h, ctx);
    fmpz_clear(omega);
}
#endif

TEST_FUNCTION_START(fq_nmod_mpoly_factor, state)
{
    slong i, j, tmul = 40;

    /* issue 2611: crash in more: branch of n_bpoly_fq_factor_lgprime */
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a;
        const char * vars[] = {"x", "y", "z"};

        fq_nmod_mpoly_ctx_init_deg(ctx, 3, ORD_LEX, 3, 1);
        fq_nmod_mpoly_init(a, ctx);

        fq_nmod_mpoly_set_str_pretty(a,
                "x^18+x^12*y^6+x^6*y^12+y^18"
                "+x^12*y^4*z^2+2*x^10*y^6*z^2+2*x^6*y^10*z^2+x^4*y^12*z^2"
                "+x^12*y^2*z^4+2*x^10*y^4*z^4+2*x^4*y^10*z^4+x^2*y^12*z^4"
                "+x^12*z^6+2*x^10*y^2*z^6+x^6*y^6*z^6+2*x^2*y^10*z^6+y^12*z^6"
                "+2*x^6*y^2*z^10+2*x^4*y^4*z^10+2*x^2*y^6*z^10"
                "+x^6*z^12+x^4*y^2*z^12+x^2*y^4*z^12+y^6*z^12+z^18",
                vars, ctx);
        check_omega(1, 1, a, ctx, fq_nmod_mpoly_factor);

        /* additional irreducible polys over GF(3) exercising the same path */
        fq_nmod_mpoly_set_str_pretty(a,
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
        check_omega(1, 1, a, ctx, fq_nmod_mpoly_factor);

        fq_nmod_mpoly_set_str_pretty(a,
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
        check_omega(1, 1, a, ctx, fq_nmod_mpoly_factor);

        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong lower;
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, t;
        slong n, nfacs, len;
        ulong expbound, powbound, pow;

        if (i % 2 == 0)
            fq_nmod_mpoly_ctx_init_rand(ctx, state, 7, 1 + n_randint(state, FLINT_BITS), 4);
        else
            fq_nmod_mpoly_ctx_init_rand(ctx, state, 2, 3, 1);

        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(t, ctx);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        nfacs = 2 + (6 + n_randint(state, 6))/n;
        powbound = 1 + n_randint(state, 3);
        powbound = 1 + n_randint(state, powbound);
        expbound = 2 + 50/nfacs/n;

        lower = 0;
        fq_nmod_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            len = 1 + n_randint(state, 1 + 20/powbound/nfacs);
            fq_nmod_mpoly_randtest_bound(t, state, len, expbound, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            pow = 1 + n_randint(state, powbound);
            if (!fq_nmod_mpoly_is_fq_nmod(t, ctx))
                lower += pow;

            fq_nmod_mpoly_pow_ui(t, t, pow, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
        }

        check_omega(lower, WORD_MAX, a, ctx, fq_nmod_mpoly_factor);

        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
