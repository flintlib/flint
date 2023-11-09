/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly.h"

/* Defined in t-gcd.c, t-gcd_brown.c, t-gcd_cofactors.c, t-gcd_hensel.c,
 * t-gcd_subresultant.c, t-gcd_zippel2.c */
#define gcd_check gcd_check_gcd_zippel2
void gcd_check(
    fmpq_mpoly_t g,
    fmpq_mpoly_t a,
    fmpq_mpoly_t b,
    const fmpq_mpoly_t gdiv,
    fmpq_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    int res;
    fmpq_mpoly_t ca, cb, cg;

    fmpq_mpoly_init(ca, ctx);
    fmpq_mpoly_init(cb, ctx);
    fmpq_mpoly_init(cg, ctx);

    res = fmpq_mpoly_gcd_zippel2(g, a, b, ctx);

    fmpq_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpq_mpoly_is_zero(gdiv, ctx))
    {
        if (!fmpq_mpoly_divides(ca, g, gdiv, ctx))
        {
            flint_printf("FAIL: Check divisor of gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    if (fmpq_mpoly_is_zero(g, ctx))
    {
        if (!fmpq_mpoly_is_zero(a, ctx) || !fmpq_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: Check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (!fmpq_mpoly_is_monic(g, ctx))
    {
        flint_printf("FAIL: Check gcd is monic\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = 1;
    res = res && fmpq_mpoly_divides(ca, a, g, ctx);
    res = res && fmpq_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        flint_printf("FAIL: Check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = fmpq_mpoly_gcd_zippel2(cg, ca, cb, ctx);
    fmpq_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd of cofactors can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpq_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: Check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fmpq_mpoly_clear(ca, ctx);
    fmpq_mpoly_clear(cb, ctx);
    fmpq_mpoly_clear(cg, ctx);
}

TEST_FUNCTION_START(fmpq_mpoly_gcd_zippel2, state)
{
    slong i, j, tmul = 10;

    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t g, a, b, t;
        const char* vars[] = {"y", "t", "x", "z"};

        fmpq_mpoly_ctx_init(ctx, 4, ORD_DEGLEX);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t, ctx);

        fmpq_mpoly_set_str_pretty(t, "x+y+z+t", vars, ctx);
        fmpq_mpoly_set_str_pretty(a, "x^2+y^2+z^2+t^2", vars, ctx);
        fmpq_mpoly_set_str_pretty(b, "x^3+y^3+z^3+t^3", vars, ctx);
        fmpq_mpoly_mul(a, a, t, ctx);
        fmpq_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 0, "example");

        fmpq_mpoly_set_str_pretty(t, "39 - t*x + 39*x^100 - t*x^101 + 39*x^3*y - t*x^4*y - 7*x^2*y^3*z^11 - 7*x^102*y^3*z^11 - 7*x^5*y^4*z^11 + 78*t^15*x^78*y^3*z^13 - 2*t^16*x^79*y^3*z^13 + x^1000*y^3*z^20 + x^1100*y^3*z^20 + x^1003*y^4*z^20 - 14*t^15*x^80*y^6*z^24 + 2*t^15*x^1078*y^6*z^33", vars, ctx);
        fmpq_mpoly_set_str_pretty(a, "39 - t*x - 7*x^2*y^3*z^11 + x^1000*y^3*z^20", vars, ctx);
        fmpq_mpoly_set_str_pretty(b, "1 + x^100 + x^3*y + 2*t^15*x^78*y^3*z^13", vars, ctx);
        fmpq_mpoly_mul(a, a, t, ctx);
        fmpq_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 0, "example");

        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, t;
        flint_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t, ctx);

        len = n_randint(state, 35) + 1;
        len1 = n_randint(state, 35) + 1;
        len2 = n_randint(state, 35) + 1;

        degbound = 2 + 100/(2*fmpq_mpoly_ctx_nvars(ctx) - 1);

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpq_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpq_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            if (fmpq_mpoly_is_zero(t, ctx))
                fmpq_mpoly_one(t, ctx);

            fmpq_mpoly_mul(a, a, t, ctx);
            fmpq_mpoly_mul(b, b, t, ctx);
            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);
            gcd_check(g, a, b, t, ctx, i, j, "sparse");
        }

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(t, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef gcd_check
