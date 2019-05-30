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
                     fmpz_mpoly_ctx_t ctx, slong i, slong j, const char * name)
{
    int res;
    fmpz_mpoly_t ca, cb, cg;

flint_printf("(%wd, %wd) %s nvars = %wd\n", i, j, name, ctx->minfo->nvars);

    fmpz_mpoly_init(ca, ctx);
    fmpz_mpoly_init(cb, ctx);
    fmpz_mpoly_init(cg, ctx);

    res = fmpz_mpoly_gcd(g, a, b, ctx);
    fmpz_mpoly_assert_canonical(g, ctx);

    if (!res)
    {
        flint_printf("Check gcd can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (fmpz_mpoly_is_zero(g, ctx))
    {
        if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check zero gcd only results from zero inputs\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
            flint_abort();
        }
        goto cleanup;
    }

    if (fmpz_sgn(g->coeffs + 0) <= 0)
    {
        printf("FAIL\n");
        flint_printf("Check gcd has positive lc\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    res = 1;
    res = res && fmpz_mpoly_divides(ca, a, g, ctx);
    res = res && fmpz_mpoly_divides(cb, b, g, ctx);
    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check divisibility\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    res = fmpz_mpoly_gcd(cg, ca, cb, ctx);
    fmpz_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        flint_abort();
    }

    if (!fmpz_mpoly_is_one(cg, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
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
    const slong max_threads = 5;
    slong i, j, k, tmul = 15;
    FLINT_TEST_INIT(state);

    flint_printf("gcd....");
    fflush(stdout);

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t g, a, b;
        const char * vars[] = {"t" ,"x", "y", "z"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(g, ctx);

        fmpz_mpoly_set_str_pretty(g, "39 - t*x + 39*x^100 - t*x^101 + 39*x^3*y - t*x^4*y - 7*x^2*y^3*z^11 - 7*x^102*y^3*z^11 - 7*x^5*y^4*z^11 + 78*t^15*x^78*y^3*z^13 - 2*t^16*x^79*y^3*z^13 + x^1000*y^3*z^20 + x^1100*y^3*z^20 + x^1003*y^4*z^20 - 14*t^15*x^80*y^6*z^24 + 2*t^15*x^1078*y^6*z^33", vars, ctx);
        fmpz_mpoly_set_str_pretty(a, "39 - t*x - 7*x^2*y^3*z^11 + x^1000*y^3*z^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "1 + x^100 + x^3*y + 2*t^15*x^78*y^3*z^13", vars, ctx);
        fmpz_mpoly_mul(a, a, g, ctx);
        fmpz_mpoly_mul(b, b, g, ctx);

        gcd_check(g, a, b, ctx, 0, 0, "example");

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    {
        int success;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t g, a, b;
        const char * vars[] = {"x" ,"y", "z", "t"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(g, ctx);

        fmpz_mpoly_set_str_pretty(a, "x^3 + 1", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "x^9999999999999999999999 + x^3333333333333333333333 + x", vars, ctx);

        flint_set_num_threads(n_randint(state, max_threads) + 1);
        success = fmpz_mpoly_gcd(g, a, b, ctx);
        if (success)
        {
            printf("FAIL\n");
            flint_printf("Check non-example\n");
            flint_abort();
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t g, a, b;
        const char * vars[] = {"x" ,"y", "z", "t"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(g, ctx);

        fmpz_mpoly_set_str_pretty(a, "(1 + x)^15*(2 + y)^18*(1 + z)^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "(2 + x)^15*(1 + y)^18*(1 - z)^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(g, "(1 - x)^15*(2 - y)^18*(1 - z)^20", vars, ctx);
        fmpz_mpoly_mul(a, a, g, ctx);
        fmpz_mpoly_mul(a, a, g, ctx);

        flint_set_num_threads(n_randint(state, max_threads) + 1);
        gcd_check(g, a, b, ctx, 0, 0, "total dense example");

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* The gcd should always work when one input is a monomial */
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

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bits(t, state, 1, coeff_bits + 1, exp_bits, ctx);
            } while (t->length != 1);
            fmpz_mpoly_randtest_bits(a, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(b, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "monomial");
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* The gcd should always work when both cofactors are monomials */
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
        len1 = n_randint(state, 25);

        exp_bits = n_randint(state, 70) + 2;
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bits(t1, state, 1, coeff_bits + 1, exp_bits1, ctx);
            } while (t1->length != 1);
            do {
                fmpz_mpoly_randtest_bits(t2, state, 1, coeff_bits + 1, exp_bits2, ctx);
            } while (t2->length != 1);
            fmpz_mpoly_randtest_bits(a, state, len1, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_mul(b, a, t1, ctx);
            fmpz_mpoly_mul(a, a, t2, ctx);

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "monomial cofactors");
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs */
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

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);

        degbound = 30/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 200);

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

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "sparse inputs");
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs with random repackings */
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

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);

        degbound = 30/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 200);

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

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "sparse input with repacking");
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs with random inflations */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t;
        mp_bitcnt_t coeff_bits;
        fmpz * shifts1, * shifts2, * strides;
        mp_bitcnt_t stride_bits, shift_bits;
        slong len, len1, len2;
        slong degbound;

        fmpz_mpoly_ctx_init_rand(ctx, state, 5);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);

        degbound = 30/(2*ctx->minfo->nvars - 1);

        coeff_bits = n_randint(state, 200);

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
            do {
                fmpz_mpoly_randtest_bound(t, state, len, coeff_bits + 1, degbound, ctx);
            } while (t->length == 0);
            fmpz_mpoly_randtest_bound(a, state, len1, coeff_bits, degbound, ctx);
            fmpz_mpoly_randtest_bound(b, state, len2, coeff_bits, degbound, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);

            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, FLINT_BITS, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts1 + k, state, shift_bits);
                fmpz_randtest_unsigned(shifts2 + k, state, shift_bits);
                fmpz_randtest_unsigned(strides + k, state, stride_bits);
            }
            fmpz_mpoly_inflate(a, a, shifts1, strides, ctx);
            fmpz_mpoly_inflate(b, b, shifts2, strides, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "sparse input with inflation");
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

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* dense inputs */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t;
        mp_bitcnt_t coeff_bits1, coeff_bits2, coeff_bits3, coeff_bits4;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        mp_bitcnt_t bits4;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        len1 = n_randint(state, 300) + 1;
        len2 = n_randint(state, 300);
        len3 = n_randint(state, 300);
        len4 = n_randint(state, 300);
 
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds1[j] = 1 + n_randint(state, 15/ctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 15/ctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 15/ctx->minfo->nvars);
        }

        bits4 = n_randint(state, 200);
        coeff_bits1 = n_randint(state, 200);
        coeff_bits2 = n_randint(state, 200);
        coeff_bits3 = n_randint(state, 200);
        coeff_bits4 = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bounds(t, state, len1, coeff_bits1 + 1, degbounds1, ctx);
            } while (t->length == 0);
            fmpz_mpoly_randtest_bounds(a, state, len2, coeff_bits2, degbounds2, ctx);
            fmpz_mpoly_randtest_bounds(b, state, len3, coeff_bits3, degbounds3, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);

            fmpz_mpoly_randtest_bits(g, state, len4, coeff_bits4, bits4, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "dense input");
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* dense inputs with repacking */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t;
        mp_limb_t rlimb;
        mp_bitcnt_t newbits;
        mp_bitcnt_t coeff_bits1, coeff_bits2, coeff_bits3, coeff_bits4;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        mp_bitcnt_t bits4;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        len1 = n_randint(state, 300) + 1;
        len2 = n_randint(state, 300);
        len3 = n_randint(state, 300);
        len4 = n_randint(state, 300);
 
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds1[j] = 1 + n_randint(state, 15/ctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 15/ctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 15/ctx->minfo->nvars);
        }

        bits4 = n_randint(state, 200);
        coeff_bits1 = n_randint(state, 200);
        coeff_bits2 = n_randint(state, 200);
        coeff_bits3 = n_randint(state, 200);
        coeff_bits4 = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            do {
                fmpz_mpoly_randtest_bounds(t, state, len1, coeff_bits1 + 1, degbounds1, ctx);
            } while (t->length == 0);
            fmpz_mpoly_randtest_bounds(a, state, len2, coeff_bits2, degbounds2, ctx);
            fmpz_mpoly_randtest_bounds(b, state, len3, coeff_bits3, degbounds3, ctx);
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

            fmpz_mpoly_randtest_bits(g, state, len4, coeff_bits4, bits4, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "dense input with repacking");
        }

        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* dense inputs with random inflations */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, g, t;
        fmpz * shifts1, * shifts2, * strides;
        mp_bitcnt_t stride_bits, shift_bits;
        mp_bitcnt_t coeff_bits1, coeff_bits2, coeff_bits3, coeff_bits4;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        mp_bitcnt_t bits4;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);

        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(t, ctx);

        len1 = n_randint(state, 300) + 1;
        len2 = n_randint(state, 300);
        len3 = n_randint(state, 300);
        len4 = n_randint(state, 300);
 
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds1[j] = 1 + n_randint(state, 15/ctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 15/ctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 15/ctx->minfo->nvars);
        }

        bits4 = n_randint(state, 200);
        coeff_bits1 = n_randint(state, 200);
        coeff_bits2 = n_randint(state, 200);
        coeff_bits3 = n_randint(state, 200);
        coeff_bits4 = n_randint(state, 200);

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
            do {
                fmpz_mpoly_randtest_bounds(t, state, len1, coeff_bits1 + 1, degbounds1, ctx);
            } while (t->length == 0);
            fmpz_mpoly_randtest_bounds(a, state, len2, coeff_bits2, degbounds2, ctx);
            fmpz_mpoly_randtest_bounds(b, state, len3, coeff_bits3, degbounds3, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
            fmpz_mpoly_mul(b, b, t, ctx);

            fmpz_mpoly_randtest_bits(g, state, len4, coeff_bits4, bits4, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts1 + k, state, shift_bits);
                fmpz_randtest_unsigned(shifts2 + k, state, shift_bits);
                fmpz_randtest_unsigned(strides + k, state, stride_bits);
            }
            fmpz_mpoly_inflate(a, a, shifts1, strides, ctx);
            fmpz_mpoly_inflate(b, b, shifts2, strides, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "dense input with inflation");
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
