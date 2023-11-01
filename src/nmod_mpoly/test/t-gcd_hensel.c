/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

/* Defined in t-gcd.c, t-gcd_brown.c, t-gcd_cofactors.c, t-gcd_hensel.c,
 * t-gcd_zippel2.c and t-gcd_zippel.c */
#define gcd_check gcd_check_gcd_hensel
void gcd_check(
    nmod_mpoly_t g,
    nmod_mpoly_t a,
    nmod_mpoly_t b,
    nmod_mpoly_t t,
    nmod_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    nmod_mpoly_t ca, cb, cg;

    nmod_mpoly_init(ca, ctx);
    nmod_mpoly_init(cb, ctx);
    nmod_mpoly_init(cg, ctx);

    if (!nmod_mpoly_gcd_hensel(g, a, b, ctx))
    {
        flint_printf("FAIL: check gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    nmod_mpoly_assert_canonical(g, ctx);

    if (nmod_mpoly_is_zero(g, ctx))
    {
        if (!nmod_mpoly_is_zero(a, ctx) || !nmod_mpoly_is_zero(b, ctx))
        {
            flint_printf("FAIL: check zero gcd\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (g->coeffs[0] != UWORD(1))
    {
        flint_printf("FAIL: check gcd is monic\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!nmod_mpoly_is_zero(t, ctx) && !nmod_mpoly_divides(cg, g, t, ctx))
    {
        flint_printf("FAIL: check gcd divisor\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!nmod_mpoly_divides(ca, a, g, ctx) ||
        !nmod_mpoly_divides(cb, b, g, ctx))
    {
        flint_printf("FAIL: check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!nmod_mpoly_gcd_hensel(cg, ca, cb, ctx))
    {
        flint_printf("FAIL: check cofactor gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    nmod_mpoly_assert_canonical(cg, ctx);

    if (!nmod_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: check gcd of cofactors is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    nmod_mpoly_clear(ca, ctx);
    nmod_mpoly_clear(cb, ctx);
    nmod_mpoly_clear(cg, ctx);
}

TEST_FUNCTION_START(nmod_mpoly_gcd_hensel, state)
{
    slong i, j;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t1, t2, t3;
        slong len, len1, len2;
        ulong degbound;
        ulong * degbounds, * degbounds1, * degbounds2;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 5, modulus);

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t1, ctx);
        nmod_mpoly_init(t2, ctx);
        nmod_mpoly_init(t3, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 50/(2*ctx->minfo->nvars - 1);
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
            nmod_mpoly_randtest_bounds(t1, state, len, degbounds, ctx);
            nmod_mpoly_randtest_bounds(t2, state, len1, degbounds1, ctx);
            nmod_mpoly_randtest_bounds(t3, state, len2, degbounds2, ctx);

            switch (n_randint(state, 4))
            {
                case 3:
                    nmod_mpoly_mul(t3, t1, t2, ctx);
                    break;
                case 2:
                    nmod_mpoly_mul(t3, t3, t1, ctx);
                    break;
                case 1:
                    nmod_mpoly_mul(t3, t3, t2, ctx);
                    break;
                default:
                    break;
            }

            nmod_mpoly_mul(a, t1, t3, ctx);
            nmod_mpoly_mul(b, t2, t3, ctx);

            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            gcd_check(g, a, b, t3, ctx, i, j, "random");
        }

        flint_free(degbounds);
        flint_free(degbounds1);
        flint_free(degbounds2);

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t1, ctx);
        nmod_mpoly_clear(t2, ctx);
        nmod_mpoly_clear(t3, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef gcd_check
