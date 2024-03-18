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
