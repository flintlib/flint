/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void gcd_check(fq_nmod_mpoly_t g, fq_nmod_mpoly_t a, fq_nmod_mpoly_t b,
                  fq_nmod_mpoly_ctx_t ctx, slong i, slong j, const char * name)
{
    int res;
    fq_nmod_mpoly_t ca, cb, cg;

    fq_nmod_mpoly_init(ca, ctx);
    fq_nmod_mpoly_init(cb, ctx);
    fq_nmod_mpoly_init(cg, ctx);

    res = fq_nmod_mpoly_gcd(g, a, b, ctx);
    fq_nmod_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("Check gcd can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (fq_nmod_mpoly_is_zero(g, ctx))
    {
        if (!fq_nmod_mpoly_is_zero(a, ctx) || !fq_nmod_mpoly_is_zero(b, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check zero gcd only results from zero inputs\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
            flint_abort();
        }
        goto cleanup;
    }

    if (!fq_nmod_is_one(g->coeffs + 0, ctx->fqctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd leading coefficient is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    res = 1;
    res = res && fq_nmod_mpoly_divides(ca, a, g, ctx);
    res = res && fq_nmod_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check divisibility\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    res = fq_nmod_mpoly_gcd(cg, ca, cb, ctx);
    fq_nmod_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (!fq_nmod_mpoly_is_one(cg, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

cleanup:

    fq_nmod_mpoly_clear(ca, ctx);
    fq_nmod_mpoly_clear(cb, ctx);
    fq_nmod_mpoly_clear(cg, ctx);
}

int
main(void)
{
    slong i, j, k, tmul = 20;
    FLINT_TEST_INIT(state);

    flint_printf("gcd....");
    fflush(stdout);

    /* The gcd should always work when one input is a monomial */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, t;
        slong len, len1, len2;
        mp_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);

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
                fq_nmod_mpoly_randtest_bits(t, state, 1, exp_bits, ctx);
            } while (t->length != 1);
            fq_nmod_mpoly_randtest_bits(a, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(b, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
            fq_nmod_mpoly_mul(b, b, t, ctx);

            fq_nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);

            gcd_check(g, a, b, ctx, i, j, "monomial");
        }

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* The gcd should always work when both cofactors are monomials */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, t1, t2;
        slong len, len1;
        mp_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t1, ctx);
        fq_nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 25);
        len1 = n_randint(state, 25);

        exp_bits = n_randint(state, 70) + 2;
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            do {
                fq_nmod_mpoly_randtest_bits(t1, state, 1, exp_bits1, ctx);
            } while (t1->length != 1);
            do {
                fq_nmod_mpoly_randtest_bits(t2, state, 1, exp_bits2, ctx);
            } while (t2->length != 1);
            fq_nmod_mpoly_randtest_bits(a, state, len1, exp_bits, ctx);
            fq_nmod_mpoly_mul(b, a, t1, ctx);
            fq_nmod_mpoly_mul(a, a, t2, ctx);

            fq_nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);

            gcd_check(g, a, b, ctx, i, j, "monomial cofactors");
        }

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t1, ctx);
        fq_nmod_mpoly_clear(t2, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }


    printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}
