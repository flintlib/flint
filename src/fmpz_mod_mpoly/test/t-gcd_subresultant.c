/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

/* Defined in t-gcd_brown.c, t-gcd_cofactors.c, t-gcd_hensel.c,
 * t-gcd_subresultant.c, t-gcd_zippel.c, t-gcd_zippel2.c */
#define gcd_check gcd_check_gcd_subresultant
void gcd_check(
    fmpz_mod_mpoly_t g,
    fmpz_mod_mpoly_t a,
    fmpz_mod_mpoly_t b,
    fmpz_mod_mpoly_t t,
    fmpz_mod_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    fmpz_mod_mpoly_t ca, cb, cg;

    fmpz_mod_mpoly_init(ca, ctx);
    fmpz_mod_mpoly_init(cb, ctx);
    fmpz_mod_mpoly_init(cg, ctx);

    if (!fmpz_mod_mpoly_gcd_subresultant(g, a, b, ctx))
    {
        flint_printf("FAIL: check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_assert_canonical(g, ctx);

    if (fmpz_mod_mpoly_is_zero(g, ctx))
    {
        if (!fmpz_mod_mpoly_is_zero(a, ctx) || !fmpz_mod_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (!fmpz_is_one(g->coeffs + 0))
    {
        flint_printf("FAIL: check gcd is monic\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mod_mpoly_is_zero(t, ctx) && !fmpz_mod_mpoly_divides(cg, g, t, ctx))
    {
        flint_printf("FAIL: check gcd divisor\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mod_mpoly_divides(ca, a, g, ctx) ||
        !fmpz_mod_mpoly_divides(cb, b, g, ctx))
    {
        flint_printf("FAIL: check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpz_mod_mpoly_gcd_subresultant(cg, ca, cb, ctx))
    {
        flint_printf("FAIL: check cofactor gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpz_mod_mpoly_assert_canonical(cg, ctx);

    if (!fmpz_mod_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fmpz_mod_mpoly_clear(ca, ctx);
    fmpz_mod_mpoly_clear(cb, ctx);
    fmpz_mod_mpoly_clear(cg, ctx);
}

TEST_FUNCTION_START(fmpz_mod_mpoly_gcd_subresultant, state)
{
    slong i, j;
    slong tmul = 10;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, t;
        slong len, len1, len2;
        slong degbound;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 3, 250);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 1 + 10/FLINT_MAX(WORD(1), ctx->minfo->nvars);

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (fmpz_mod_mpoly_is_zero(t, ctx))
                fmpz_mod_mpoly_one(t, ctx);
            fmpz_mod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            fmpz_mod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            fmpz_mod_mpoly_mul(a, a, t, ctx);
            fmpz_mod_mpoly_mul(b, b, t, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            gcd_check(g, a, b, t, ctx, i, j, "random small");
        }

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef gcd_check
