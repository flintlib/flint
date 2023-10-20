/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly_factor.h"

/* Defined t-factor.c, t-factor_wang.c, t-factor_zassenhaus.c,
 * t-factor_zippel.c */
#ifndef check_omega
#define check_omega check_omega
void check_omega(
        slong lower,
        slong upper,
        const fmpz_mpoly_t p,
        const fmpz_mpoly_ctx_t ctx,
        int factor_fun(fmpz_mpoly_factor_t, const fmpz_mpoly_t, const fmpz_mpoly_ctx_t))
{
    slong i;
    fmpz_mpoly_t q;
    fmpz_mpoly_factor_t g, h;
    fmpz_t omega;

    fmpz_init(omega);
    fmpz_mpoly_factor_init(g, ctx);
    fmpz_mpoly_factor_init(h, ctx);
    fmpz_mpoly_init(q, ctx);

    if (!factor_fun(g, p, ctx))
    {
        flint_printf("check factorization 1 could be computed\n");
        fflush(stdout);
        flint_abort();
    }

    if (factor_fun != fmpz_mpoly_factor)
    {
        if (!fmpz_mpoly_factor(h, p, ctx))
        {
            flint_printf("check factorization 2 could be computed\n");
            fflush(stdout);
            flint_abort();
        }
    }

    for (i = 0; i < g->num; i++)
    {
        if (g->poly[i].length < 1 || fmpz_sgn(g->poly[i].coeffs + 0) <= 0)
        {
            flint_printf("factorization is not unit normal\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, lower) < 0 || fmpz_cmp_si(omega, upper) > 0)
    {
        flint_printf("factorization has wrong number of factors\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_mpoly_factor_expand(q, g, ctx);

    if (!fmpz_mpoly_equal(q, p, ctx))
    {
        flint_printf("factorization does not match original polynomial\n");
        fflush(stdout);
        flint_abort();
    }

    if (factor_fun != fmpz_mpoly_factor)
    {
        fmpz_mpoly_factor_sort(g, ctx);
        fmpz_mpoly_factor_sort(h, ctx);
        if (fmpz_mpoly_factor_cmp(g, h, ctx) != 0)
        {
            flint_printf("factorizations do not match\n");
            fflush(stdout);
            flint_abort();
        }
    }

    for (i = 0; i < g->num; i++)
    {
        fmpz_mpoly_factor(h, g->poly + i, ctx);
        if (h->num != 1 || !fmpz_is_one(h->exp + 0))
        {
            flint_printf("FAIL:\nfactor is reducible\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_mpoly_clear(q, ctx);
    fmpz_mpoly_factor_clear(g, ctx);
    fmpz_mpoly_factor_clear(h, ctx);
    fmpz_clear(omega);
    return;
}
#endif

TEST_FUNCTION_START(fmpz_mpoly_factor_zippel, state)
{
    slong i, j, tmul = 20;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong lower;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, t;
        flint_bitcnt_t coeff_bits;
        slong n, nfacs, len;
        ulong expbound, powbound, pow;

        fmpz_mpoly_ctx_init_rand(ctx, state, 7);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(t, ctx);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        nfacs = 1 + (5 + n_randint(state, 5))/n;
        expbound = 3 + 40/nfacs/n;
        powbound = 1 + n_randint(state, 3);

        lower = 0;
        fmpz_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            do {
                len = 1 + n_randint(state, 7);
                coeff_bits = 10 + n_randint(state, 1000)/nfacs;
                fmpz_mpoly_randtest_bound(t, state, len, coeff_bits, expbound, ctx);
            } while (t->length == 0);
            pow = 1 + n_randint(state, powbound);
            if (!fmpz_mpoly_is_fmpz(t, ctx))
                lower += pow;
            fmpz_mpoly_pow_ui(t, t, pow, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
        }

        check_omega(lower, WORD_MAX, a, ctx, fmpz_mpoly_factor_zippel);

        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
