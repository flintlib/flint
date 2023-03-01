/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

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

    if (!fmpz_mod_mpoly_gcd_hensel(g, a, b, ctx))
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

    if (!fmpz_mod_mpoly_gcd_hensel(cg, ca, cb, ctx))
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


int
main(void)
{
    slong i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_hensel....");
    fflush(stdout);

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t a, b, g, t1, t2, t3;
        slong len, len1, len2;
        ulong degbound;
        ulong * degbounds, * degbounds1, * degbounds2;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 5, 100);

        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(a, ctx);
        fmpz_mod_mpoly_init(b, ctx);
        fmpz_mod_mpoly_init(t1, ctx);
        fmpz_mod_mpoly_init(t2, ctx);
        fmpz_mod_mpoly_init(t3, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 42/(2*ctx->minfo->nvars - 1);
        degbounds = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        degbounds1 = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));
        degbounds2 = (ulong *) flint_malloc(ctx->minfo->nvars*sizeof(ulong));

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds[j] = n_randint(state, degbound + 1) + 1;
            degbounds1[j] = n_randint(state, degbound + 1) + 1;
            degbounds2[j] = n_randint(state, degbound + 1) + 1;
        }

        for (j = 0; j < 6; j++)
        {
            fmpz_mod_mpoly_randtest_bounds(t1, state, len, degbounds, ctx);
            fmpz_mod_mpoly_randtest_bounds(t2, state, len1, degbounds1, ctx);
            fmpz_mod_mpoly_randtest_bounds(t3, state, len2, degbounds2, ctx);

            switch (n_randint(state, 4))
            {
                case 3:
                    fmpz_mod_mpoly_mul(t3, t1, t2, ctx);
                    break;
                case 2:
                    fmpz_mod_mpoly_mul(t3, t3, t1, ctx);
                    break;
                case 1:
                    fmpz_mod_mpoly_mul(t3, t3, t2, ctx);
                    break;
                default:
                    break;
            }

            fmpz_mod_mpoly_mul(a, t1, t3, ctx);
            fmpz_mod_mpoly_mul(b, t2, t3, ctx);

            fmpz_mod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            gcd_check(g, a, b, t3, ctx, i, j, "random");
        }

        flint_free(degbounds);
        flint_free(degbounds1);
        flint_free(degbounds2);

        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(a, ctx);
        fmpz_mod_mpoly_clear(b, ctx);
        fmpz_mod_mpoly_clear(t1, ctx);
        fmpz_mod_mpoly_clear(t2, ctx);
        fmpz_mod_mpoly_clear(t3, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    flint_printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
