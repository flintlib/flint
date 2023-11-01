/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly_factor.h"

void check_omega(slong lower, slong upper, const fmpq_mpoly_t p, const fmpq_mpoly_ctx_t ctx)
{
    slong i;
    fmpq_t t;
    fmpq_mpoly_t q;
    fmpq_mpoly_factor_t g, h;
    fmpz_t omega;

    fmpq_init(t);
    fmpz_init(omega);
    fmpq_mpoly_factor_init(g, ctx);
    fmpq_mpoly_factor_init(h, ctx);
    fmpq_mpoly_init(q, ctx);

    if (!fmpq_mpoly_factor(g, p, ctx))
    {
        flint_printf("check factorization could be computed\n");
        fflush(stdout);
        flint_abort();
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

    fmpq_mpoly_factor_expand(q, g, ctx);
    if (!fmpq_mpoly_equal(q, p, ctx))
    {
        flint_printf("factorization does not match original polynomial\n");
        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < g->num; i++)
    {
        fmpq_mpoly_factor(h, g->poly + i, ctx);
        if (h->num != 1 || !fmpz_is_one(h->exp + 0))
        {
            flint_printf("FAIL:\nfactor is reducible\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpq_mpoly_factor_make_monic(g, ctx);
    for (i = 0; i < g->num; i++)
    {
        if (fmpq_mpoly_length(g->poly + i, ctx) < 1 ||
            (fmpq_mpoly_get_term_coeff_fmpq(t, g->poly + i, 0, ctx),
                                                             !fmpq_is_one(t)))
        {
            flint_printf("monic factorization is not monic\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpq_mpoly_factor_expand(q, g, ctx);
    if (!fmpq_mpoly_equal(q, p, ctx))
    {
        flint_printf("monic factorization does not match original polynomial\n");
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_factor_make_integral(g, ctx);
    for (i = 0; i < g->num; i++)
    {
        if (fmpq_mpoly_length(g->poly + i, ctx) < 1 ||
            (fmpq_mpoly_content(t, g->poly + i, ctx), !fmpq_is_one(t)) ||
            (fmpq_mpoly_get_term_coeff_fmpq(t, g->poly + i, 0, ctx),
                                                            fmpq_sgn(t) <= 0))
        {
            flint_printf("integral factorization is not integral\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpq_mpoly_factor_expand(q, g, ctx);
    if (!fmpq_mpoly_equal(q, p, ctx))
    {
        flint_printf("integral factorization does not match original polynomial\n");
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_clear(q, ctx);
    fmpq_mpoly_factor_clear(g, ctx);
    fmpq_mpoly_factor_clear(h, ctx);
    fmpz_clear(omega);
    fmpq_clear(t);
}

TEST_FUNCTION_START(fmpq_mpoly_factor, state)
{
    slong i, j, tmul = 25;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong lower;
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, t;
        flint_bitcnt_t coeff_bits;
        slong n, nfacs, len;
        ulong expbound, powbound, pow;

        fmpq_mpoly_ctx_init_rand(ctx, state, 8);

        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(t, ctx);

        n = FLINT_MAX(WORD(1), ctx->zctx->minfo->nvars);
        nfacs = 1 + (5 + n_randint(state, 5))/n;
        expbound = 3 + 40/nfacs/n;
        powbound = 1 + n_randint(state, 3);

        lower = 0;
        fmpq_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            do {
                len = 1 + n_randint(state, 7);
                coeff_bits = 10 + n_randint(state, 200)/nfacs;
                fmpq_mpoly_randtest_bound(t, state, len, coeff_bits, expbound, ctx);
            } while (fmpq_mpoly_length(t, ctx) < 1);
            pow = 1 + n_randint(state, powbound);
            if (!fmpq_mpoly_is_fmpq(t, ctx))
                lower += pow;
            fmpq_mpoly_pow_ui(t, t, pow, ctx);
            fmpq_mpoly_mul(a, a, t, ctx);
        }

        check_omega(lower, WORD_MAX, a, ctx);

        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
