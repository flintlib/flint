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

TEST_FUNCTION_START(fmpz_mpoly_factor, state)
{
    slong i, j, tmul = 25;

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;

        fmpz_mpoly_ctx_init(ctx, 8, ORD_LEX);
        fmpz_mpoly_init(f, ctx);

        fmpz_mpoly_set_str_pretty(f,
                "x1^809*x2^75*x3^384*x4^324*x5^74*x6^788*x7^83*x8^414+"
                "x1^805*x2^343*x3^595*x4^246*x5^32*x6^90*x7^473*x8^591+"
                "x1^718*x2^108*x3^680*x4^368*x5^358*x8^276+"
                "x1^683*x2^533*x4^649*x5^619*x6^136*x7^223*x8^610+"
                "x2^617*x3^777*x4^799*x5^443*x6^545*x7^166*x8^216+"
                "x1^485*x2^646*x3^424*x4^265*x5^416*x6^400*x7^278+"
                "x1^336*x2^149*x3^361*x4^691*x5^629*x6^282*x7^530*x8^259+"
                "x1^266*x3^258*x5^422*x6^637*x7^244*x8^236+"
                "x1^74*x2^812*x3^162*x4^417*x5^71*x6^188*x7^258*x8^637+"
                "x1^37*x2^604*x3^94*x4^474*x6^853*x7^521*x8^250", NULL, ctx);

        check_omega(1, 1, f, ctx, fmpz_mpoly_factor);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong lower;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, t;
        flint_bitcnt_t coeff_bits;
        slong n, nfacs, len;
        ulong expbound, powbound, pow;

        fmpz_mpoly_ctx_init_rand(ctx, state, 8);

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

        check_omega(lower, WORD_MAX, a, ctx, fmpz_mpoly_factor);

        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
