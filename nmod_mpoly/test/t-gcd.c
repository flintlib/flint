/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void gcd_check(nmod_mpoly_t g, nmod_mpoly_t a, nmod_mpoly_t b,
                     nmod_mpoly_ctx_t ctx, slong i, slong j, const char * name)
{
    int res;
    nmod_mpoly_t ca, cb, cg;

    nmod_mpoly_init(ca, ctx);
    nmod_mpoly_init(cb, ctx);
    nmod_mpoly_init(cg, ctx);

    res = nmod_mpoly_gcd(g, a, b, ctx);
    nmod_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("Check gcd can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (nmod_mpoly_is_zero(g, ctx))
    {
        if (!nmod_mpoly_is_zero(a, ctx) || !nmod_mpoly_is_zero(b, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check zero gcd only results from zero inputs\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
            flint_abort();
        }
        goto cleanup;
    }

    if (g->coeffs[0] != UWORD(1))
    {
        printf("FAIL\n");
        flint_printf("Check gcd leading coefficient is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    res = 1;
    res = res && nmod_mpoly_divides(ca, a, g, ctx);
    res = res && nmod_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check divisibility\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    res = nmod_mpoly_gcd(cg, ca, cb, ctx);
    nmod_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (!nmod_mpoly_is_one(cg, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

cleanup:

    nmod_mpoly_clear(ca, ctx);
    nmod_mpoly_clear(cb, ctx);
    nmod_mpoly_clear(cg, ctx);
}

int
main(void)
{
    slong i, j, tmul = 20;

    FLINT_TEST_INIT(state);

    flint_printf("gcd....");
    fflush(stdout);

    /* The gcd should always work when one input is a monomial */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t;
        slong len, len1, len2;
        mp_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);

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
                nmod_mpoly_randtest_bits(t, state, 1, exp_bits, ctx);
            } while (t->length != 1);
            nmod_mpoly_randtest_bits(a, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(b, state, len2, exp_bits2, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_mul(b, b, t, ctx);

            nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);

            gcd_check(g, a, b, ctx, i, j, "monomial");
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }


    for (i = 0; i < 2 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, ca, cb, cg, t;
        slong len, len1, len2;
        slong degbound;
        mp_limb_t modulus;
        int res;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, modulus < 3000 ? 4 : 5, modulus);

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(ca, ctx);
        nmod_mpoly_init(cb, ctx);
        nmod_mpoly_init(cg, ctx);
        nmod_mpoly_init(t, ctx);

        len = n_randint(state, 100) + 1;
        len1 = n_randint(state, 200);
        len2 = n_randint(state, 200);

        degbound = 100/ctx->minfo->nvars/ctx->minfo->nvars;

        for (j = 0; j < 4; j++)
        {
            do {
                nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            } while (t->length == 0);
            nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_mul(b, b, t, ctx);

            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            res = nmod_mpoly_gcd(g, a, b, ctx);
            if (!res) {
                continue;
            }
            nmod_mpoly_assert_canonical(g, ctx);

            if (nmod_mpoly_is_zero(g, ctx))
            {
                if (!nmod_mpoly_is_zero(a, ctx) || !nmod_mpoly_is_zero(b, ctx))
                {
                    printf("FAIL\n");
                    flint_printf("Check zero gcd only results from zero inputs\ni = %wd, j = %wd\n", i ,j);
                    flint_abort();
                }
                continue;
            }

            if (g->coeffs[0] != UWORD(1))
            {
                printf("FAIL\n");
                flint_printf("Check gcd is monic\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = 1;
            res = res && nmod_mpoly_divides(ca, a, g, ctx);
            res = res && nmod_mpoly_divides(cb, b, g, ctx);
            if (!res)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            res = nmod_mpoly_gcd(cg, ca, cb, ctx);

            if (!res)
                continue;

            if (!nmod_mpoly_equal_ui(cg, UWORD(1), ctx))
            {
                printf("FAIL\n");
                flint_printf("Check cofactors are relatively prime\ni = %wd, j = %wd\n", i ,j);                
                flint_abort();
            }
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(ca, ctx);
        nmod_mpoly_clear(cb, ctx);
        nmod_mpoly_clear(cg, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }


    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

