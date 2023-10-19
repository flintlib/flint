/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpq_mpoly.h"

/* Defined in t-gcd.c, t-gcd_brown.c, t-gcd_cofactors.c, t-gcd_hensel.c,
 * t-gcd_subresultant.c, t-gcd_zippel2.c */
#define gcd_check gcd_check_gcd_cofactors
void gcd_check(
    fmpq_mpoly_t g,
    fmpq_mpoly_t abar,
    fmpq_mpoly_t bbar,
    fmpq_mpoly_t a,
    fmpq_mpoly_t b,
    const fmpq_mpoly_t gdiv,
    fmpq_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    int res;
    fmpq_mpoly_t ca, cb, cg, u, v, w;

    fmpq_mpoly_init(ca, ctx);
    fmpq_mpoly_init(cb, ctx);
    fmpq_mpoly_init(cg, ctx);
    fmpq_mpoly_init(u, ctx);
    fmpq_mpoly_init(v, ctx);
    fmpq_mpoly_init(w, ctx);

    res = fmpq_mpoly_gcd_cofactors(g, abar, bbar, a, b, ctx);

    fmpq_mpoly_assert_canonical(g, ctx);
    fmpq_mpoly_assert_canonical(abar, ctx);
    fmpq_mpoly_assert_canonical(bbar, ctx);

    if (!res)
    {
        flint_printf("Check gcd can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpq_mpoly_is_zero(gdiv, ctx))
    {
        if (!fmpq_mpoly_divides(ca, g, gdiv, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check divisor of gcd\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    fmpq_mpoly_mul(ca, g, abar, ctx);
    fmpq_mpoly_mul(cb, g, bbar, ctx);
    if (!fmpq_mpoly_equal(ca, a, ctx) || !fmpq_mpoly_equal(cb, b, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check cofactors i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (fmpq_mpoly_is_zero(g, ctx))
    {
        if (!fmpq_mpoly_is_zero(a, ctx) || !fmpq_mpoly_is_zero(b, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check zero gcd only results from zero inputs\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (!fmpq_mpoly_is_monic(g, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd is monic\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(u, b, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, a, u, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (a, u): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(v, b, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, a, v, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (a, v): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(w, b, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, a, w, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (a, w): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(u, a, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, u, b, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (u, b): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(v, a, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, v, b, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (v, b): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(w, a, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, w, b, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (w, b): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(u, a, ctx);
    fmpq_mpoly_set(v, b, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, u, v, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (u, v): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(v, a, ctx);
    fmpq_mpoly_set(u, b, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, v, u, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (v, u): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(u, a, ctx);
    fmpq_mpoly_set(w, b, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, u, w, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (u, w): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(w, a, ctx);
    fmpq_mpoly_set(u, b, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, w, u, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (w, u): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(v, a, ctx);
    fmpq_mpoly_set(w, b, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, v, w, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (v, w): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_set(w, a, ctx);
    fmpq_mpoly_set(v, b, ctx);
    fmpq_mpoly_gcd_cofactors(u, v, w, w, v, ctx);
    fmpq_mpoly_assert_canonical(u, ctx);
    fmpq_mpoly_assert_canonical(v, ctx);
    fmpq_mpoly_assert_canonical(w, ctx);
    if (!fmpq_mpoly_equal(g, u, ctx) || !fmpq_mpoly_equal(abar, v, ctx) || !fmpq_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (w, v): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = fmpq_mpoly_gcd_cofactors(cg, ca, cb, abar, bbar, ctx);
    fmpq_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpq_mpoly_is_one(cg, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpq_mpoly_equal(ca, abar, ctx) || !fmpq_mpoly_equal(cb, bbar, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check cofactors of cofactors\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = fmpq_mpoly_gcd_cofactors(cg, abar, bbar, abar, bbar, ctx);
    fmpq_mpoly_assert_canonical(cg, ctx);

    if (!fmpq_mpoly_equal(ca, abar, ctx) || !fmpq_mpoly_equal(cb, bbar, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check cofactors of cofactors with aliasing\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fmpq_mpoly_clear(ca, ctx);
    fmpq_mpoly_clear(cb, ctx);
    fmpq_mpoly_clear(cg, ctx);
    fmpq_mpoly_clear(u, ctx);
    fmpq_mpoly_clear(v, ctx);
    fmpq_mpoly_clear(w, ctx);
}

TEST_FUNCTION_START(fmpq_mpoly_gcd_cofactors, state)
{
    const slong max_threads = 5;
    slong i, j, k, tmul = 3;

    {
        int success;
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t g, abar, bbar, a, b;
        const char * vars[] = {"x" ,"y", "z", "t"};

        fmpq_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(abar, ctx);
        fmpq_mpoly_init(bbar, ctx);

        fmpq_mpoly_set_str_pretty(a, "x^3 + 1", vars, ctx);
        fmpq_mpoly_set_str_pretty(b, "x^9999999999999999999999 + x^3333333333333333333333 + x", vars, ctx);

        flint_set_num_threads(n_randint(state, max_threads) + 1);
        success = fmpq_mpoly_gcd_cofactors(g, abar, bbar, a, b, ctx);
        if (success)
        {
            printf("FAIL\n");
            flint_printf("Check non-example\n");
            fflush(stdout);
            flint_abort();
        }

        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(abar, ctx);
        fmpq_mpoly_clear(bbar, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* The gcd should always work when one input is a monomial */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, abar, bbar, t;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(abar, ctx);
        fmpq_mpoly_init(bbar, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t, ctx);

        len = n_randint(state, 25);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);
        if (n_randlimb(state) & UWORD(1))
            len1 = FLINT_MIN(len1, WORD(1));
        else
            len2 = FLINT_MIN(len2, WORD(1));

        exp_bits = n_randint(state, 70) + 2;
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpq_mpoly_randtest_bits(t, state, 1, coeff_bits + 1, exp_bits, ctx);
            } while (t->zpoly->length != 1);
            fmpq_mpoly_randtest_bits(a, state, len1, coeff_bits, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(b, state, len2, coeff_bits, exp_bits2, ctx);
            fmpq_mpoly_mul(a, a, t, ctx);
            fmpq_mpoly_mul(b, b, t, ctx);

            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "monomial");
        }

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(abar, ctx);
        fmpq_mpoly_clear(bbar, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* The gcd should always work when both cofactors are monomials */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, abar, bbar, t1, t2;
        slong len, len1;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(abar, ctx);
        fmpq_mpoly_init(bbar, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t1, ctx);
        fmpq_mpoly_init(t2, ctx);

        len = n_randint(state, 25);
        len1 = n_randint(state, 25);

        exp_bits = n_randint(state, 70) + 2;
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpq_mpoly_randtest_bits(t1, state, 1, coeff_bits + 1, exp_bits1, ctx);
            } while (t1->zpoly->length != 1);
            do {
                fmpq_mpoly_randtest_bits(t2, state, 1, coeff_bits + 1, exp_bits2, ctx);
            } while (t2->zpoly->length != 1);
            fmpq_mpoly_randtest_bits(a, state, len1, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_mul(b, a, t1, ctx);
            fmpq_mpoly_mul(t2, a, t2, ctx);

            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, t2, b, a, ctx, i, j, "monomial cofactors");
        }

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(abar, ctx);
        fmpq_mpoly_clear(bbar, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(t1, ctx);
        fmpq_mpoly_clear(t2, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* one input divides the other */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_t c;
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, abar, bbar, t1, t2;
        slong len, len1, len2;
        mp_limb_t exp_bound, exp_bound1, exp_bound2;
        flint_bitcnt_t coeff_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_init(c);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(abar, ctx);
        fmpq_mpoly_init(bbar, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t1, ctx);
        fmpq_mpoly_init(t2, ctx);

        len = n_randint(state, 5);
        len1 = n_randint(state, 5);
        len2 = n_randint(state, 5);

        exp_bound = n_randint(state, 100) + 2;
        exp_bound1 = n_randint(state, 100) + 2;
        exp_bound2 = n_randint(state, 100) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(t1, state, len1, coeff_bits + 1, exp_bound1, ctx);
            fmpq_mpoly_randtest_bound(t2, state, len2, coeff_bits + 1, exp_bound2, ctx);
            fmpq_mpoly_mul(b, t1, t2, ctx);
            fmpz_randtest(c, state, coeff_bits + 1);
            fmpq_mpoly_scalar_mul_fmpz(a, t2, c, ctx);
            fmpz_randtest(c, state, coeff_bits + 1);
            fmpq_mpoly_scalar_mul_fmpz(b, b, c, ctx);

            fmpq_mpoly_randtest_bound(g, state, len, coeff_bits, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            if ((j%2) == 0)
                fmpq_mpoly_swap(a, b, ctx);

            gcd_check(g, abar, bbar, a, b, t2, ctx, i, j, "one input divides the other");
        }

        fmpz_clear(c);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(abar, ctx);
        fmpq_mpoly_clear(bbar, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(t1, ctx);
        fmpq_mpoly_clear(t2, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, abar, bbar, t;
        flint_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;

        fmpq_mpoly_ctx_init_rand(ctx, state, 5);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(abar, ctx);
        fmpq_mpoly_init(bbar, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t, ctx);

        len = n_randint(state, 10) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 25/(2*ctx->zctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpq_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            } while (t->zpoly->length == 0);
            fmpq_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpq_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpq_mpoly_mul(a, a, t, ctx);
            fmpq_mpoly_mul(b, b, t, ctx);

            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "sparse inputs");
        }

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(abar, ctx);
        fmpq_mpoly_clear(bbar, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs with random repackings */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, abar, bbar, t;
        mp_limb_t rlimb;
        flint_bitcnt_t coeff_bits, newbits;
        slong len, len1, len2;
        slong degbound;

        fmpq_mpoly_ctx_init_rand(ctx, state, 5);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(abar, ctx);
        fmpq_mpoly_init(bbar, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t, ctx);

        len = n_randint(state, 10) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 25/(2*ctx->zctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpq_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            } while (t->zpoly->length == 0);
            fmpq_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpq_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpq_mpoly_mul(a, a, t, ctx);
            fmpq_mpoly_mul(b, b, t, ctx);

            rlimb = n_randlimb(state);

            if (rlimb & UWORD(3))
            {
                newbits = a->zpoly->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->zctx->minfo);
                fmpq_mpoly_repack_bits(a, a, newbits, ctx);
            }

            if (rlimb & UWORD(12))
            {
                newbits = b->zpoly->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->zctx->minfo);
                fmpq_mpoly_repack_bits(b, b, newbits, ctx);
            }

            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "sparse input with repacking");
        }

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(abar, ctx);
        fmpq_mpoly_clear(bbar, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs with random inflations */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, abar, bbar, t;
        flint_bitcnt_t coeff_bits;
        fmpz * shifts1, * shifts2, * strides;
        flint_bitcnt_t stride_bits, shift_bits;
        slong len, len1, len2;
        slong degbound;

        fmpq_mpoly_ctx_init_rand(ctx, state, 5);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(abar, ctx);
        fmpq_mpoly_init(bbar, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t, ctx);

        len = n_randint(state, 10) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 25/(2*ctx->zctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 200);

        stride_bits = n_randint(state, 100) + 2;
        shift_bits = n_randint(state, 100) + 2;

        shifts1 = flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz));
        shifts2 = flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz));
        strides = flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz));
        for (k = 0; k < ctx->zctx->minfo->nvars; k++)
        {
            fmpz_init(shifts1 + k);
            fmpz_init(shifts2 + k);
            fmpz_init(strides + k);
        }

        for (j = 0; j < 4; j++)
        {
            do {
                fmpq_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            } while (t->zpoly->length == 0);
            fmpq_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpq_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpq_mpoly_mul(a, a, t, ctx);
            fmpq_mpoly_mul(b, b, t, ctx);

            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            for (k = 0; k < ctx->zctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts1 + k, state, shift_bits);
                fmpz_randtest_unsigned(shifts2 + k, state, shift_bits);
                fmpz_randtest_unsigned(strides + k, state, stride_bits);
            }
            fmpq_mpoly_inflate(a, a, shifts1, strides, ctx);
            fmpq_mpoly_inflate(b, b, shifts2, strides, ctx);

            for (k = 0; k < ctx->zctx->minfo->nvars; k++)
            {
                if (fmpz_cmp(shifts1 + k, shifts2 + k) > 0)
                    fmpz_set(shifts1 + k, shifts2 + k);
            }
            fmpq_mpoly_inflate(t, t, shifts1, strides, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "sparse input with inflation");
        }

        for (k = 0; k < ctx->zctx->minfo->nvars; k++)
        {
            fmpz_clear(shifts1 + k);
            fmpz_clear(shifts2 + k);
            fmpz_clear(strides + k);
        }
        flint_free(shifts1);
        flint_free(shifts2);
        flint_free(strides);

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(abar, ctx);
        fmpq_mpoly_clear(bbar, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* dense inputs */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, abar, bbar, t;
        flint_bitcnt_t coeff_bits1, coeff_bits2, coeff_bits3, coeff_bits4;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;

        fmpq_mpoly_ctx_init_rand(ctx, state, 4);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(abar, ctx);
        fmpq_mpoly_init(bbar, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t, ctx);

        len1 = n_randint(state, 150) + 1;
        len2 = n_randint(state, 150);
        len3 = n_randint(state, 150);
        len4 = n_randint(state, 150);

        for (j = 0; j < ctx->zctx->minfo->nvars; j++)
        {
            degbounds1[j] = 1 + n_randint(state, 12/ctx->zctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 12/ctx->zctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 12/ctx->zctx->minfo->nvars);
        }

        bits4 = n_randint(state, 100);
        coeff_bits1 = n_randint(state, 100);
        coeff_bits2 = n_randint(state, 100);
        coeff_bits3 = n_randint(state, 100);
        coeff_bits4 = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpq_mpoly_randtest_bounds(t, state, len1, coeff_bits1 + 1, degbounds1, ctx);
            } while (t->zpoly->length == 0);
            fmpq_mpoly_randtest_bounds(a, state, len2, coeff_bits2, degbounds2, ctx);
            fmpq_mpoly_randtest_bounds(b, state, len3, coeff_bits3, degbounds3, ctx);
            fmpq_mpoly_mul(a, a, t, ctx);
            fmpq_mpoly_mul(b, b, t, ctx);

            fmpq_mpoly_randtest_bits(g, state, len4, coeff_bits4, bits4, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "dense input");
        }

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(abar, ctx);
        fmpq_mpoly_clear(bbar, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* dense inputs with repacking */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, abar, bbar, t;
        mp_limb_t rlimb;
        flint_bitcnt_t newbits;
        flint_bitcnt_t coeff_bits1, coeff_bits2, coeff_bits3, coeff_bits4;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;

        fmpq_mpoly_ctx_init_rand(ctx, state, 4);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(abar, ctx);
        fmpq_mpoly_init(bbar, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t, ctx);

        len1 = n_randint(state, 150) + 1;
        len2 = n_randint(state, 150);
        len3 = n_randint(state, 150);
        len4 = n_randint(state, 150);

        for (j = 0; j < ctx->zctx->minfo->nvars; j++)
        {
            degbounds1[j] = 1 + n_randint(state, 12/ctx->zctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 12/ctx->zctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 12/ctx->zctx->minfo->nvars);
        }

        bits4 = n_randint(state, 100);
        coeff_bits1 = n_randint(state, 100);
        coeff_bits2 = n_randint(state, 100);
        coeff_bits3 = n_randint(state, 100);
        coeff_bits4 = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpq_mpoly_randtest_bounds(t, state, len1, coeff_bits1 + 1, degbounds1, ctx);
            } while (t->zpoly->length == 0);
            fmpq_mpoly_randtest_bounds(a, state, len2, coeff_bits2, degbounds2, ctx);
            fmpq_mpoly_randtest_bounds(b, state, len3, coeff_bits3, degbounds3, ctx);
            fmpq_mpoly_mul(a, a, t, ctx);
            fmpq_mpoly_mul(b, b, t, ctx);

            rlimb = n_randlimb(state);

            if (rlimb & UWORD(3))
            {
                newbits = a->zpoly->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->zctx->minfo);
                fmpq_mpoly_repack_bits(a, a, newbits, ctx);
            }

            if (rlimb & UWORD(12))
            {
                newbits = b->zpoly->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->zctx->minfo);
                fmpq_mpoly_repack_bits(b, b, newbits, ctx);
            }

            fmpq_mpoly_randtest_bits(g, state, len4, coeff_bits4, bits4, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "dense input with repacking");
        }

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(abar, ctx);
        fmpq_mpoly_clear(bbar, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* dense inputs with random inflations */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, abar, bbar, t;
        fmpz * shifts1, * shifts2, * strides;
        flint_bitcnt_t stride_bits, shift_bits;
        flint_bitcnt_t coeff_bits1, coeff_bits2, coeff_bits3, coeff_bits4;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;

        fmpq_mpoly_ctx_init_rand(ctx, state, 4);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(abar, ctx);
        fmpq_mpoly_init(bbar, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t, ctx);

        len1 = n_randint(state, 150) + 1;
        len2 = n_randint(state, 150);
        len3 = n_randint(state, 150);
        len4 = n_randint(state, 150);

        for (j = 0; j < ctx->zctx->minfo->nvars; j++)
        {
            degbounds1[j] = 1 + n_randint(state, 12/ctx->zctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 12/ctx->zctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 12/ctx->zctx->minfo->nvars);
        }

        bits4 = n_randint(state, 100);
        coeff_bits1 = n_randint(state, 100);
        coeff_bits2 = n_randint(state, 100);
        coeff_bits3 = n_randint(state, 100);
        coeff_bits4 = n_randint(state, 100);

        stride_bits = n_randint(state, 100) + 2;
        shift_bits = n_randint(state, 100) + 2;

        shifts1 = flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz));
        shifts2 = flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz));
        strides = flint_malloc(ctx->zctx->minfo->nvars*sizeof(fmpz));
        for (k = 0; k < ctx->zctx->minfo->nvars; k++)
        {
            fmpz_init(shifts1 + k);
            fmpz_init(shifts2 + k);
            fmpz_init(strides + k);
        }

        for (j = 0; j < 4; j++)
        {
            do {
                fmpq_mpoly_randtest_bounds(t, state, len1, coeff_bits1 + 1, degbounds1, ctx);
            } while (t->zpoly->length == 0);
            fmpq_mpoly_randtest_bounds(a, state, len2, coeff_bits2, degbounds2, ctx);
            fmpq_mpoly_randtest_bounds(b, state, len3, coeff_bits3, degbounds3, ctx);
            fmpq_mpoly_mul(a, a, t, ctx);
            fmpq_mpoly_mul(b, b, t, ctx);

            fmpq_mpoly_randtest_bits(g, state, len4, coeff_bits4, bits4, ctx);

            for (k = 0; k < ctx->zctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts1 + k, state, shift_bits);
                fmpz_randtest_unsigned(shifts2 + k, state, shift_bits);
                fmpz_randtest_unsigned(strides + k, state, stride_bits);
            }
            fmpq_mpoly_inflate(a, a, shifts1, strides, ctx);
            fmpq_mpoly_inflate(b, b, shifts2, strides, ctx);

            for (k = 0; k < ctx->zctx->minfo->nvars; k++)
            {
                if (fmpz_cmp(shifts1 + k, shifts2 + k) > 0)
                    fmpz_set(shifts1 + k, shifts2 + k);
            }
            fmpq_mpoly_inflate(t, t, shifts1, strides, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "dense input with inflation");
        }

        for (k = 0; k < ctx->zctx->minfo->nvars; k++)
        {
            fmpz_clear(shifts1 + k);
            fmpz_clear(shifts2 + k);
            fmpz_clear(strides + k);
        }
        flint_free(shifts1);
        flint_free(shifts2);
        flint_free(strides);

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(abar, ctx);
        fmpq_mpoly_clear(bbar, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef gcd_check
