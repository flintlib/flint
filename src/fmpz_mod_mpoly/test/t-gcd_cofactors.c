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
#include "fmpz_mod_mpoly.h"

/* Defined in t-gcd_brown.c, t-gcd_cofactors.c, t-gcd_hensel.c,
 * t-gcd_subresultant.c, t-gcd_zippel.c, t-gcd_zippel2.c */
#define gcd_check gcd_check_gcd_cofactors
void gcd_check(
    fmpz_mod_mpoly_t g,
    fmpz_mod_mpoly_t abar,
    fmpz_mod_mpoly_t bbar,
    fmpz_mod_mpoly_t a,
    fmpz_mod_mpoly_t b,
    const fmpz_mod_mpoly_t gdiv,
    fmpz_mod_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    int res;
    fmpz_mod_mpoly_t ca, cb, cg, u, v, w;

    fmpz_mod_mpoly_init(ca, ctx);
    fmpz_mod_mpoly_init(cb, ctx);
    fmpz_mod_mpoly_init(cg, ctx);
    fmpz_mod_mpoly_init(u, ctx);
    fmpz_mod_mpoly_init(v, ctx);
    fmpz_mod_mpoly_init(w, ctx);

    res = fmpz_mod_mpoly_gcd_cofactors(g, abar, bbar, a, b, ctx);

    fmpz_mod_mpoly_assert_canonical(g, ctx);
    fmpz_mod_mpoly_assert_canonical(abar, ctx);
    fmpz_mod_mpoly_assert_canonical(bbar, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mod_mpoly_is_zero(gdiv, ctx))
    {
        if (!fmpz_mod_mpoly_divides(ca, g, gdiv, ctx))
        {
            flint_printf("FAIL: Check divisor of gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_mod_mpoly_mul(ca, g, abar, ctx);
    fmpz_mod_mpoly_mul(cb, g, bbar, ctx);
    if (!fmpz_mod_mpoly_equal(ca, a, ctx) || !fmpz_mod_mpoly_equal(cb, b, ctx))
    {
        flint_printf("FAIL: Check cofactors\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (fmpz_mod_mpoly_is_zero(g, ctx))
    {
        if (!fmpz_mod_mpoly_is_zero(a, ctx) || !fmpz_mod_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: Check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (!fmpz_is_one(g->coeffs + 0))
    {
        flint_printf("FAIL: Check gcd is monic\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(u, b, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, a, u, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (a, u): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(v, b, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, a, v, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (a, v): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(w, b, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, a, w, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (a, w): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(u, a, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, u, b, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (u, b): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(v, a, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, v, b, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (v, b): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(w, a, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, w, b, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (w, b): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(u, a, ctx);
    fmpz_mod_mpoly_set(v, b, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, u, v, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (u, v): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(v, a, ctx);
    fmpz_mod_mpoly_set(u, b, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, v, u, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (v, u): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(u, a, ctx);
    fmpz_mod_mpoly_set(w, b, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, u, w, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (u, w): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(w, a, ctx);
    fmpz_mod_mpoly_set(u, b, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, w, u, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (w, u): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(v, a, ctx);
    fmpz_mod_mpoly_set(w, b, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, v, w, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (v, w): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_set(w, a, ctx);
    fmpz_mod_mpoly_set(v, b, ctx);
    fmpz_mod_mpoly_gcd_cofactors(u, v, w, w, v, ctx);
    if (!fmpz_mod_mpoly_equal(g, u, ctx) || !fmpz_mod_mpoly_equal(abar, v, ctx) || !fmpz_mod_mpoly_equal(bbar, w, ctx))
    {
        flint_printf("FAIL (u, v, w), (w, v): i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = fmpz_mod_mpoly_gcd_cofactors(cg, ca, cb, abar, bbar, ctx);
    fmpz_mod_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd of cofactors can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mod_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: Check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mod_mpoly_equal(ca, abar, ctx) || !fmpz_mod_mpoly_equal(cb, bbar, ctx))
    {
        flint_printf("FAIL: Check cofactors of cofactors\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = fmpz_mod_mpoly_gcd_cofactors(cg, abar, bbar, abar, bbar, ctx);
    fmpz_mod_mpoly_assert_canonical(cg, ctx);

    if (!fmpz_mod_mpoly_equal(ca, abar, ctx) || !fmpz_mod_mpoly_equal(cb, bbar, ctx))
    {
        flint_printf("FAIL: Check cofactors of cofactors with aliasing\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fmpz_mod_mpoly_clear(ca, ctx);
    fmpz_mod_mpoly_clear(cb, ctx);
    fmpz_mod_mpoly_clear(cg, ctx);
    fmpz_mod_mpoly_clear(u, ctx);
    fmpz_mod_mpoly_clear(v, ctx);
    fmpz_mod_mpoly_clear(w, ctx);
}

TEST_FUNCTION_START(fmpz_mod_mpoly_gcd_cofactors, state)
{
    const slong max_threads = 5;
    slong i, j, k, tmul = 2;
    fmpz_t p;

    fmpz_init_set_ui(p, 1);
    fmpz_mul_2exp(p, p, 100);
    fmpz_nextprime(p, p, 1);

    for (i = 3; i <= 6; i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t g, abar, bbar, a, b, t;

        fmpz_mod_mpoly_ctx_init(ctx, i, ORD_DEGREVLEX, p);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        fmpz_mod_mpoly_one(g, ctx);
        fmpz_mod_mpoly_one(a, ctx);
        fmpz_mod_mpoly_one(b, ctx);
        for (j = 0; j < i; j++)
        {
            fmpz_mod_mpoly_gen(t, j, ctx);
            fmpz_mod_mpoly_add_si(t, t, 1, ctx);
            fmpz_mod_mpoly_mul(g, g, t, ctx);
            fmpz_mod_mpoly_gen(t, j, ctx);
            fmpz_mod_mpoly_sub_si(t, t, 2, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_gen(t, j, ctx);
            fmpz_mod_mpoly_add_si(t, t, 2, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);
        }
        fmpz_mod_mpoly_sub_si(g, g, 2, ctx);
        fmpz_mod_mpoly_add_si(a, a, 2, ctx);
        fmpz_mod_mpoly_sub_si(b, b, 2, ctx);

        fmpz_mod_mpoly_mul(a, a, g, ctx);
        fmpz_mod_mpoly_mul(b, b, g, ctx);
        fmpz_mod_mpoly_set(t, g, ctx);

        gcd_check(g, abar, bbar, a, b, t, ctx, i, 0, "dense examples");

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    if (0) {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t g, abar, bbar, a, b, t;
        const char * vars[] = {"t" ,"x", "y", "z"};

        fmpz_mod_mpoly_ctx_init(ctx, 4, ORD_LEX, p);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        fmpz_mod_mpoly_set_str_pretty(t, "39 - t*x + 39*x^100 - t*x^101 + 39*x^3*y - t*x^4*y - 7*x^2*y^3*z^11 - 7*x^102*y^3*z^11 - 7*x^5*y^4*z^11 + 78*t^15*x^78*y^3*z^13 - 2*t^16*x^79*y^3*z^13 + x^1000*y^3*z^20 + x^1100*y^3*z^20 + x^1003*y^4*z^20 - 14*t^15*x^80*y^6*z^24 + 2*t^15*x^1078*y^6*z^33", vars, ctx);
        fmpz_mod_mpoly_set_str_pretty(a, "39 - t*x - 7*x^2*y^3*z^11 + x^1000*y^3*z^20", vars, ctx);
        fmpz_mod_mpoly_set_str_pretty(b, "1 + x^100 + x^3*y + 2*t^15*x^78*y^3*z^13", vars, ctx);
        fmpz_mod_mpoly_mul(a, a, t, ctx);
        fmpz_mod_mpoly_mul(b, b, t, ctx);

        gcd_check(g, abar, bbar, a, b, t, ctx, 0, 0, "example");

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    {
        int success;
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t g, abar, bbar, a, b;
        const char * vars[] = {"x" ,"y", "z", "t"};

        fmpz_mod_mpoly_ctx_init(ctx, 4, ORD_LEX, p);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);

        fmpz_mod_mpoly_set_str_pretty(a, "x^3 + 1", vars, ctx);
        fmpz_mod_mpoly_set_str_pretty(b, "x^9999999999999999999999 + x^3333333333333333333333 + x", vars, ctx);

        flint_set_num_threads(n_randint(state, max_threads) + 1);
        success = fmpz_mod_mpoly_gcd_cofactors(g, abar, bbar, a, b, ctx);
        if (success)
        {
            flint_printf("FAIL\n");
            flint_printf("Check non-example\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t g, abar, bbar, a, b, t;
        const char * vars[] = {"x" ,"y", "z", "t"};

        fmpz_mod_mpoly_ctx_init(ctx, 4, ORD_LEX, p);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        fmpz_mod_mpoly_set_str_pretty(a, "(1 + x)^1*(2 + y)^1*(1 + z)^2", vars, ctx);
        fmpz_mod_mpoly_set_str_pretty(b, "(2 + x)^1*(1 + y)^1*(1 - z)^2", vars, ctx);
        fmpz_mod_mpoly_set_str_pretty(t, "(1 - x)^1*(2 - y)^1*(1 - z)^2", vars, ctx);
        fmpz_mod_mpoly_mul(a, a, t, ctx);
        fmpz_mod_mpoly_mul(b, b, t, ctx);

        flint_set_num_threads(n_randint(state, max_threads) + 1);
        gcd_check(g, abar, bbar, a, b, t, ctx, 0, 0, "total dense example");

        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* The gcd should always work when one input is a monomial */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, abar, bbar, t;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 10, 150);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);

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

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mod_mpoly_randtest_bits(t, state, 1, exp_bits, ctx);
            } while (t->length != 1);
            fmpz_mod_mpoly_randtest_bits(a, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bits(b, state, len2, exp_bits2, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);

            fmpz_mod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "monomial");
        }

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* The gcd should always work when both cofactors are monomials */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, abar, bbar, t1, t2;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 10, 150);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t1, ctx);
        fmpz_mod_mpoly_init(t2, ctx);

        len = n_randint(state, 25);
        len1 = n_randint(state, 25);

        exp_bits = n_randint(state, 70) + 2;
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mod_mpoly_randtest_bits(t1, state, 1, exp_bits1, ctx);
            } while (t1->length != 1);
            do {
                fmpz_mod_mpoly_randtest_bits(t2, state, 1, exp_bits2, ctx);
            } while (t2->length != 1);
            fmpz_mod_mpoly_randtest_bits(a, state, len1, exp_bits, ctx);
            fmpz_mod_mpoly_mul(b, a, t1, ctx);
            fmpz_mod_mpoly_mul(t2, a, t2, ctx);

            fmpz_mod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, t2, b, a, ctx, i, j, "monomial cofactors");
        }

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t1, ctx);
        fmpz_mod_mpoly_clear(t2, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* one input divides the other */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_t c;
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, abar, bbar, t1, t2;
        slong len, len1, len2;
        mp_limb_t exp_bound, exp_bound1, exp_bound2;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 10, 150);

        fmpz_init(c);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t1, ctx);
        fmpz_mod_mpoly_init(t2, ctx);

        len = n_randint(state, 5);
        len1 = n_randint(state, 5);
        len2 = n_randint(state, 5);

        exp_bound = n_randint(state, 100) + 2;
        exp_bound1 = n_randint(state, 100) + 2;
        exp_bound2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(t1, state, len1, exp_bound1, ctx);
            fmpz_mod_mpoly_randtest_bound(t2, state, len2, exp_bound2, ctx);
            fmpz_mod_mpoly_mul(b, t1, t2, ctx);
            fmpz_randm(c, state, fmpz_mod_mpoly_ctx_modulus(ctx));
            fmpz_mod_mpoly_scalar_mul_fmpz(a, t2, c, ctx);
            fmpz_randm(c, state, fmpz_mod_mpoly_ctx_modulus(ctx));
            fmpz_mod_mpoly_scalar_mul_fmpz(b, b, c, ctx);

            fmpz_mod_mpoly_randtest_bound(g, state, len, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            if ((j%2) == 0)
                fmpz_mod_mpoly_swap(a, b, ctx);

            gcd_check(g, abar, bbar, a, b, t2, ctx, i, j, "one input divides the other");
        }

        fmpz_clear(c);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t1, ctx);
        fmpz_mod_mpoly_clear(t2, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, abar, bbar, t;
        slong len, len1, len2;
        slong degbound;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 5, 150);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 25/(2*ctx->minfo->nvars - 1);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            } while (t->length == 0);
            fmpz_mod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            fmpz_mod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);

            fmpz_mod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "sparse inputs");
        }

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs with random repackings */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, abar, bbar, t;
        mp_limb_t rlimb;
        flint_bitcnt_t newbits;
        slong len, len1, len2;
        slong degbound;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 5, 150);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 25/(2*ctx->minfo->nvars - 1);

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (fmpz_mod_mpoly_is_zero(t, ctx))
                fmpz_mod_mpoly_one(t, ctx);
            fmpz_mod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            fmpz_mod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);

            rlimb = n_randlimb(state);

            if (rlimb & UWORD(3))
            {
                newbits = a->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                fmpz_mod_mpoly_repack_bits(a, a, newbits, ctx);
            }

            if (rlimb & UWORD(12))
            {
                newbits = b->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                fmpz_mod_mpoly_repack_bits(b, b, newbits, ctx);
            }

            fmpz_mod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "sparse input with repacking");
        }

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs with random inflations */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, abar, bbar, t;
        fmpz * shifts1, * shifts2, * strides;
        flint_bitcnt_t stride_bits, shift_bits;
        slong len, len1, len2;
        slong degbound;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 5, 150);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 25/(2*ctx->minfo->nvars - 1);

        stride_bits = n_randint(state, 100) + 2;
        shift_bits = n_randint(state, 100) + 2;

        shifts1 = flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        shifts2 = flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        strides = flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_init(shifts1 + k);
            fmpz_init(shifts2 + k);
            fmpz_init(strides + k);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (fmpz_mod_mpoly_is_one(t, ctx))
                fmpz_mod_mpoly_one(t, ctx);
            fmpz_mod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            fmpz_mod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);

            fmpz_mod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts1 + k, state, shift_bits);
                fmpz_randtest_unsigned(shifts2 + k, state, shift_bits);
                fmpz_randtest_unsigned(strides + k, state, stride_bits);
            }
            fmpz_mod_mpoly_inflate(a, a, shifts1, strides, ctx);
            fmpz_mod_mpoly_inflate(b, b, shifts2, strides, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                if (fmpz_cmp(shifts1 + k, shifts2 + k) > 0)
                    fmpz_set(shifts1 + k, shifts2 + k);
            }
            fmpz_mod_mpoly_inflate(t, t, shifts1, strides, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "sparse input with inflation");
        }

        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_clear(shifts1 + k);
            fmpz_clear(shifts2 + k);
            fmpz_clear(strides + k);
        }
        flint_free(shifts1);
        flint_free(shifts2);
        flint_free(strides);

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* dense inputs */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, abar, bbar, t;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 4, 150);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        len1 = n_randint(state, 150) + 1;
        len2 = n_randint(state, 150);
        len3 = n_randint(state, 150);
        len4 = n_randint(state, 150);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds1[j] = 1 + n_randint(state, 12/ctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 12/ctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 12/ctx->minfo->nvars);
        }

        bits4 = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bounds(t, state, len1, degbounds1, ctx);
            if (fmpz_mod_mpoly_is_zero(t, ctx))
                fmpz_mod_mpoly_one(t, ctx);
            fmpz_mod_mpoly_randtest_bounds(a, state, len2, degbounds2, ctx);
            fmpz_mod_mpoly_randtest_bounds(b, state, len3, degbounds3, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);

            fmpz_mod_mpoly_randtest_bits(g, state, len4, bits4, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "dense input");
        }

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* dense inputs with repacking */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, abar, bbar, t;
        mp_limb_t rlimb;
        flint_bitcnt_t newbits;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 4, 150);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        len1 = n_randint(state, 150) + 1;
        len2 = n_randint(state, 150);
        len3 = n_randint(state, 150);
        len4 = n_randint(state, 150);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds1[j] = 1 + n_randint(state, 12/ctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 12/ctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 12/ctx->minfo->nvars);
        }

        bits4 = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bounds(t, state, len1, degbounds1, ctx);
            if (fmpz_mod_mpoly_is_zero(t, ctx))
                fmpz_mod_mpoly_one(t, ctx);
            fmpz_mod_mpoly_randtest_bounds(a, state, len2, degbounds2, ctx);
            fmpz_mod_mpoly_randtest_bounds(b, state, len3, degbounds3, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);

            rlimb = n_randlimb(state);

            if (rlimb & UWORD(3))
            {
                newbits = a->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                fmpz_mod_mpoly_repack_bits(a, a, newbits, ctx);
            }

            if (rlimb & UWORD(12))
            {
                newbits = b->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                fmpz_mod_mpoly_repack_bits(b, b, newbits, ctx);
            }

            fmpz_mod_mpoly_randtest_bits(g, state, len4, bits4, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "dense input with repacking");
        }

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* dense inputs with random inflations */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, abar, bbar, t;
        fmpz * shifts1, * shifts2, * strides;
        flint_bitcnt_t stride_bits, shift_bits;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 4, 150);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(abar, ctx);
        fmpz_mod_mpoly_init(bbar, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        len1 = n_randint(state, 150) + 1;
        len2 = n_randint(state, 150);
        len3 = n_randint(state, 150);
        len4 = n_randint(state, 150);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds1[j] = 1 + n_randint(state, 12/ctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 12/ctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 12/ctx->minfo->nvars);
        }

        bits4 = n_randint(state, 100);

        stride_bits = n_randint(state, 100) + 2;
        shift_bits = n_randint(state, 100) + 2;

        shifts1 = flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        shifts2 = flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        strides = flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_init(shifts1 + k);
            fmpz_init(shifts2 + k);
            fmpz_init(strides + k);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bounds(t, state, len1, degbounds1, ctx);
            if (fmpz_mod_mpoly_is_zero(t, ctx))
                fmpz_mod_mpoly_one(t, ctx);
            fmpz_mod_mpoly_randtest_bounds(a, state, len2, degbounds2, ctx);
            fmpz_mod_mpoly_randtest_bounds(b, state, len3, degbounds3, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);

            fmpz_mod_mpoly_randtest_bits(g, state, len4, bits4, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts1 + k, state, shift_bits);
                fmpz_randtest_unsigned(shifts2 + k, state, shift_bits);
                fmpz_randtest_unsigned(strides + k, state, stride_bits);
            }
            fmpz_mod_mpoly_inflate(a, a, shifts1, strides, ctx);
            fmpz_mod_mpoly_inflate(b, b, shifts2, strides, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                if (fmpz_cmp(shifts1 + k, shifts2 + k) > 0)
                    fmpz_set(shifts1 + k, shifts2 + k);
            }
            fmpz_mod_mpoly_inflate(t, t, shifts1, strides, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, abar, bbar, a, b, t, ctx, i, j, "dense input with inflation");
        }

        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_clear(shifts1 + k);
            fmpz_clear(shifts2 + k);
            fmpz_clear(strides + k);
        }
        flint_free(shifts1);
        flint_free(shifts2);
        flint_free(strides);

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(abar, ctx);
        fmpz_mod_mpoly_clear(bbar, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    fmpz_clear(p);

    TEST_FUNCTION_END(state);
}
#undef gcd_check
