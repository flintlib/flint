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
#define gcd_check gcd_check_gcd_hensel
void gcd_check(
    fmpq_mpoly_t g,
    fmpq_mpoly_t a,
    fmpq_mpoly_t b,
    fmpq_mpoly_t t,
    fmpq_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    fmpq_mpoly_t ca, cb, cg;

    fmpq_mpoly_init(ca, ctx);
    fmpq_mpoly_init(cb, ctx);
    fmpq_mpoly_init(cg, ctx);

    if (!fmpq_mpoly_gcd_hensel(g, a, b, ctx))
    {
        flint_printf("FAIL: check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_assert_canonical(g, ctx);

    if (fmpq_mpoly_is_zero(g, ctx))
    {
        if (!fmpq_mpoly_is_zero(a, ctx) || !fmpq_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (!fmpq_mpoly_is_monic(g, ctx))
    {
        flint_printf("FAIL: check gcd is monic\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpq_mpoly_is_zero(t, ctx) && !fmpq_mpoly_divides(cg, g, t, ctx))
    {
        flint_printf("FAIL: check gcd divisor\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpq_mpoly_divides(ca, a, g, ctx) ||
        !fmpq_mpoly_divides(cb, b, g, ctx))
    {
        flint_printf("FAIL: check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fmpq_mpoly_gcd_hensel(cg, ca, cb, ctx))
    {
        flint_printf("FAIL: check cofactor gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    fmpq_mpoly_assert_canonical(cg, ctx);

    if (!fmpq_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fmpq_mpoly_clear(ca, ctx);
    fmpq_mpoly_clear(cb, ctx);
    fmpq_mpoly_clear(cg, ctx);
}

TEST_FUNCTION_START(fmpq_mpoly_gcd_hensel, state)
{
    slong i, j, tmul = 10;

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t a, b, g, t1, t2, t3;
        slong len, len1, len2;
        ulong degbound;
        ulong * degbounds, * degbounds1, * degbounds2;
        flint_bitcnt_t coeff_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 5);

        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(a, ctx);
        fmpq_mpoly_init(b, ctx);
        fmpq_mpoly_init(t1, ctx);
        fmpq_mpoly_init(t2, ctx);
        fmpq_mpoly_init(t3, ctx);

        degbound = 42/(2*fmpq_mpoly_ctx_nvars(ctx) - 1);
        degbounds = FLINT_ARRAY_ALLOC(fmpq_mpoly_ctx_nvars(ctx), ulong);
        degbounds1 = FLINT_ARRAY_ALLOC(fmpq_mpoly_ctx_nvars(ctx), ulong);
        degbounds2 = FLINT_ARRAY_ALLOC(fmpq_mpoly_ctx_nvars(ctx), ulong);

        for (j = 0; j < fmpq_mpoly_ctx_nvars(ctx); j++)
        {
            degbounds[j] = n_randint(state, degbound + 1) + 1;
            degbounds1[j] = n_randint(state, degbound + 1) + 1;
            degbounds2[j] = n_randint(state, degbound + 1) + 1;
        }

        for (j = 0; j < 6; j++)
        {
            len = n_randint(state, 10) + 1;
            len1 = n_randint(state, 10);
            len2 = n_randint(state, 15);

            coeff_bits = n_randint(state, 200) + 10;
            fmpq_mpoly_randtest_bounds(t1, state, coeff_bits, len, degbounds, ctx);
            coeff_bits = n_randint(state, 200) + 10;
            fmpq_mpoly_randtest_bounds(t2, state, coeff_bits, len1, degbounds1, ctx);
            coeff_bits = n_randint(state, 200) + 10;
            fmpq_mpoly_randtest_bounds(t3, state, coeff_bits, len2, degbounds2, ctx);

            switch (n_randint(state, 4))
            {
                case 3:
                    fmpq_mpoly_mul(t3, t1, t2, ctx);
                    break;
                case 2:
                    fmpq_mpoly_mul(t3, t3, t1, ctx);
                    break;
                case 1:
                    fmpq_mpoly_mul(t3, t3, t2, ctx);
                    break;
                default:
                    break;
            }

            fmpq_mpoly_mul(a, t1, t3, ctx);
            fmpq_mpoly_mul(b, t2, t3, ctx);

            coeff_bits = n_randint(state, 300) + 10;
            fmpq_mpoly_randtest_bits(g, state, coeff_bits, len, FLINT_BITS, ctx);

            if (fmpq_mpoly_length(a, ctx) < 4000 && fmpq_mpoly_length(b, ctx) < 4000)
                gcd_check(g, a, b, t3, ctx, i, j, "random");
        }

        flint_free(degbounds);
        flint_free(degbounds1);
        flint_free(degbounds2);

        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(a, ctx);
        fmpq_mpoly_clear(b, ctx);
        fmpq_mpoly_clear(t1, ctx);
        fmpq_mpoly_clear(t2, ctx);
        fmpq_mpoly_clear(t3, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef gcd_check
