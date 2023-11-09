/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

/* Defined in t-gcd.c, t-gcd_brown.c, t-gcd_cofactors.c, t-gcd_hensel.c,
 * t-gcd_subresultant.c, t-gcd_zippel2.c */
#define gcd_check gcd_check_gcd_zippel2
void gcd_check(
    fmpz_mpoly_t g,
    fmpz_mpoly_t a,
    fmpz_mpoly_t b,
    const fmpz_mpoly_t gdiv,
    fmpz_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    int res;
    fmpz_mpoly_t ca, cb, cg;

    fmpz_mpoly_init(ca, ctx);
    fmpz_mpoly_init(cb, ctx);
    fmpz_mpoly_init(cg, ctx);

    res = fmpz_mpoly_gcd_zippel2(g, a, b, ctx);

    fmpz_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mpoly_is_zero(gdiv, ctx))
    {
        if (!fmpz_mpoly_divides(ca, g, gdiv, ctx))
        {
            flint_printf("FAIL: Check divisor of gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    if (fmpz_mpoly_is_zero(g, ctx))
    {
        if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: Check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (fmpz_sgn(g->coeffs + 0) <= 0)
    {
        flint_printf("FAIL: Check gcd has positive lc\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = 1;
    res = res && fmpz_mpoly_divides(ca, a, g, ctx);
    res = res && fmpz_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        flint_printf("FAIL: Check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    res = fmpz_mpoly_gcd_zippel2(cg, ca, cb, ctx);
    fmpz_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        flint_printf("FAIL: Check gcd of cofactors can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: Check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fmpz_mpoly_clear(ca, ctx);
    fmpz_mpoly_clear(cb, ctx);
    fmpz_mpoly_clear(cg, ctx);
}

TEST_FUNCTION_START(fmpz_mpoly_gcd_zippel2, state)
{
    slong i, j, tmul = 20;

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t g, a, b, t;
        const char* vars[] = {"y", "t", "x", "z"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_DEGLEX);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        fmpz_mpoly_set_str_pretty(t, "x+y+z+t", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "x^2+y^2+z^2+t^2", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "x^3+y^3+z^3+t^3", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 0, "example");

        fmpz_mpoly_set_str_pretty(t, "39 - t*x + 39*x^100 - t*x^101 + 39*x^3*y - t*x^4*y - 7*x^2*y^3*z^11 - 7*x^102*y^3*z^11 - 7*x^5*y^4*z^11 + 78*t^15*x^78*y^3*z^13 - 2*t^16*x^79*y^3*z^13 + x^1000*y^3*z^20 + x^1100*y^3*z^20 + x^1003*y^4*z^20 - 14*t^15*x^80*y^6*z^24 + 2*t^15*x^1078*y^6*z^33", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "39 - t*x - 7*x^2*y^3*z^11 + x^1000*y^3*z^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "1 + x^100 + x^3*y + 2*t^15*x^78*y^3*z^13", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 0, "example");

        fmpz_mpoly_set_str_pretty(t, "y + t^2 + x^3 + z^4", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "y*t + 1 + (x - z^5)*(y + t)", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "y*t + 1 + (x - z^5)*(y - t + x)", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 1, "trigger unlucky ksub");

        fmpz_mpoly_set_str_pretty(t, "y + 33857*t^2 + 35153*x^3 + 40433*z^4", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "y^4 + t^3 + x^2 + 1", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "y^3 + t^4 + x^2 + z", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 2, "trigger zipple no match");

        fmpz_mpoly_set_str_pretty(t, "y + t + x^3 + z^3", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "(x - z^4)*y + 1", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "(x + z)*y + t", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 3, "trigger ksub lc kill");

        fmpz_mpoly_set_str_pretty(t, "y + t + x^3 + z^3", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "(x - z^4 + 33857*(x*z))*y + 1", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "(x + z)*y + t", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 4, "trigger ksub lc kill mod p");

        fmpz_mpoly_set_str_pretty(t, "(1 + x^10)*t*y + t + x + z", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "(1 + x + t)*(x - 33857*x^2 + 35153*z^4)*y*t + z + x*y", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "(2*x - t^2)*(x - 33857*x^2 + 35153*z^4)*y^2*t + t*z + y", vars, ctx);
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 5, "trigger gcd lc terms vanish mod p");

        if (FLINT_BITS == 64)
        {
            fmpz_mpoly_set_str_pretty(t, "t*y + t + x^9999999999 + z^9999999999", vars, ctx);
            fmpz_mpoly_set_str_pretty(a, "t + y + x^9999999999", vars, ctx);
            fmpz_mpoly_set_str_pretty(b, "t^2 + t*z^9999999999 + y + 1", vars, ctx);
        }
        else
        {
            fmpz_mpoly_set_str_pretty(t, "t*y + t + x^999999 + z^9999999", vars, ctx);
            fmpz_mpoly_set_str_pretty(a, "t + y + x^999999", vars, ctx);
            fmpz_mpoly_set_str_pretty(b, "t^2 + t*z^999999 + y + 1", vars, ctx);
        }
        fmpz_mpoly_mul(a, a, t, ctx);
        fmpz_mpoly_mul(b, b, t, ctx);
        gcd_check(g, a, b, t, ctx, -1, 6, "trigger big p");

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t;
        flint_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        if (ctx->minfo->nvars < 3)
        {
            fmpz_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 35) + 1;
        len1 = n_randint(state, 35) + 1;
        len2 = n_randint(state, 35) + 1;

        degbound = 2 + 100/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            if (fmpz_mpoly_is_zero(t, ctx))
                fmpz_mpoly_one(t, ctx);

            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);
            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);
            gcd_check(g, a, b, t, ctx, i, j, "sparse");
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef gcd_check
