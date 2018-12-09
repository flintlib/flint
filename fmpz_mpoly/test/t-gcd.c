/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

void gcd_check(fmpz_mpoly_t g, fmpz_mpoly_t a, fmpz_mpoly_t b,
                               fmpz_mpoly_ctx_t ctx, slong i, slong j, slong k)
{
    int res;
    fmpz_mpoly_t ca, cb, cg;

    fmpz_mpoly_init(ca, ctx);
    fmpz_mpoly_init(cb, ctx);
    fmpz_mpoly_init(cg, ctx);

    res = fmpz_mpoly_gcd(g, a, b, ctx);
    fmpz_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("Check gcd can be computed\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
        flint_abort();
    }

    if (fmpz_mpoly_is_zero(g, ctx))
    {
        if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check zero gcd only results from zero inputs\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
            flint_abort();
        }
        goto cleanup;
    }

    if (fmpz_sgn(g->coeffs + 0) <= 0)
    {
        printf("FAIL\n");
        flint_printf("Check gcd has positive lc\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
        flint_abort();
    }

    res = 1;
    res = res && fmpz_mpoly_divides_monagan_pearce(ca, a, g, ctx);
    res = res && fmpz_mpoly_divides_monagan_pearce(cb, b, g, ctx);
    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check divisibility\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
        flint_abort();
    }

    res = fmpz_mpoly_gcd(cg, ca, cb, ctx);
    fmpz_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors can be computed\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
        flint_abort();
    }

    if (!fmpz_mpoly_is_one(cg, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors is one\n"
                                       "i = %wd, j = %wd, k = %wd\n", i, j, k);
        flint_abort();
    }

cleanup:

    fmpz_mpoly_clear(ca, ctx);
    fmpz_mpoly_clear(cb, ctx);
    fmpz_mpoly_clear(cg, ctx);
}

int
main(void)
{
    int i, j, tmul = 30;
    FLINT_TEST_INIT(state);

    flint_printf("gcd....");
    fflush(stdout);

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t g, a, b;
        const char * vars[] = {"x" ,"y", "z", "t"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(g, ctx);

        fmpz_mpoly_set_str_pretty(g, "39 - t*x + 39*x^100 - t*x^101 + 39*x^3*y - t*x^4*y - 7*x^2*y^3*z^11 - 7*x^102*y^3*z^11 - 7*x^5*y^4*z^11 + 78*t^15*x^78*y^3*z^13 - 2*t^16*x^79*y^3*z^13 + x^1000*y^3*z^20 + x^1100*y^3*z^20 + x^1003*y^4*z^20 - 14*t^15*x^80*y^6*z^24 + 2*t^15*x^1078*y^6*z^33", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "39 - t*x - 7*x^2*y^3*z^11 + x^1000*y^3*z^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "1 + x^100 + x^3*y + 2*t^15*x^78*y^3*z^13", vars, ctx);
        fmpz_mpoly_mul(a, a, g, ctx);
        fmpz_mpoly_mul(b, b, g, ctx);

        gcd_check(g, a, b, ctx, 0,0,0);

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* k = 1: The gcd should always work when one input is a monomial */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t;
        slong len, len1, len2;
        mp_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

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

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bound(t, state, 1, coeff_bits + 1, exp_bits, ctx);
            } while (t->length != 1);
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);

            gcd_check(g, a, b, ctx, i, j, 1);
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* k = 2: The gcd should always work when both cofactors are monomials */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t1, t2;
        slong len, len1;
        mp_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 25);
        len1 = n_randint(state, 5);

        exp_bits = n_randint(state, 70) + 2;
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bits(t1, state, 1, coeff_bits + 1, exp_bits1, ctx);
            } while (t1->length != 1);
            do {
                fmpz_mpoly_randtest_bits(t2, state, 1, coeff_bits + 1, exp_bits2, ctx);
            } while (t2->length != 1);
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_mul(b, a, t1, ctx);
            fmpz_mpoly_mul(a, a, t2, ctx);

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);

            gcd_check(g, a, b, ctx, i, j, 2);
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* k = 3: test dense inputs */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t;
        mp_bitcnt_t coeff_bits;
        slong len, len1, len2;
        slong degbound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 5);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 25) + 1;
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        degbound = 25/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 40);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            } while (t->length == 0);
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            gcd_check(g, a, b, ctx, i, j, 3);
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* k = 4: test dense inputs with random repackings */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t;
        mp_limb_t rlimb;
        mp_bitcnt_t coeff_bits, newbits;
        slong len, len1, len2;
        slong degbound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 5);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 25) + 1;
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        degbound = 25/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 40);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            } while (t->length == 0);
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);

            rlimb = n_randlimb(state);

            if (rlimb & UWORD(3))
            {
                newbits = a->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                fmpz_mpoly_repack_bits(a, a, newbits, ctx);
            }

            if (rlimb & UWORD(12))
            {
                newbits = b->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                fmpz_mpoly_repack_bits(b, b, newbits, ctx);
            }

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            gcd_check(g, a, b, ctx, i, j, 4);
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

