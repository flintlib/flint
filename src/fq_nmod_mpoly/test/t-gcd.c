/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

/* Defined in t-gcd_brown.c, t-gcd.c, t-gcd_cofactors.c, t-gcd_hensel.c,
 * t-gcd_zippel2.c, t-gcd_zippel.c */
#define gcd_check gcd_check_gcd
void gcd_check(
    fq_nmod_mpoly_t g,
    fq_nmod_mpoly_t a,
    fq_nmod_mpoly_t b,
    const fq_nmod_mpoly_t gdiv,
    const fq_nmod_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
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
        fflush(stdout);
        flint_abort();
    }

    if (!fq_nmod_mpoly_is_zero(gdiv, ctx))
    {
        if (!fq_nmod_mpoly_divides(ca, g, gdiv, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check divisor of gcd\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    if (fq_nmod_mpoly_is_zero(g, ctx))
    {
        if (!fq_nmod_mpoly_is_zero(a, ctx) || !fq_nmod_mpoly_is_zero(b, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check zero gcd only results from zero inputs\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
        goto cleanup;
    }

    if (!fq_nmod_mpoly_is_monic(g, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd leading coefficient is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
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
        fflush(stdout);
        flint_abort();
    }

    res = fq_nmod_mpoly_gcd(cg, ca, cb, ctx);
    fq_nmod_mpoly_assert_canonical(cg, ctx);

    if (!res)
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors can be computed\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!fq_nmod_mpoly_is_one(cg, ctx))
    {
        printf("FAIL\n");
        flint_printf("Check gcd of cofactors is one\n"
                                         "i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    fq_nmod_mpoly_clear(ca, ctx);
    fq_nmod_mpoly_clear(cb, ctx);
    fq_nmod_mpoly_clear(cg, ctx);
}

TEST_FUNCTION_START(fq_nmod_mpoly_gcd, state)
{
    slong i, j, k, tmul = 7;

    /* The gcd should always work when one input is a monomial */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, t;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

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
            fq_nmod_mpoly_randtest_bits(t, state, 1, exp_bits, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            fq_nmod_mpoly_randtest_bits(a, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(b, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
            fq_nmod_mpoly_mul(b, b, t, ctx);

            fq_nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);

            gcd_check(g, a, b, t, ctx, i, j, "monomial");
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
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

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
            fq_nmod_mpoly_randtest_bits(t1, state, 1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(t2, state, 1, exp_bits2, ctx);
            if (t1->length != 1 || t2->length != 1)
            {
                flint_printf("FAIL:\ncheck random monomial generation\n");
                fflush(stdout);
                flint_abort();
            }
            fq_nmod_mpoly_randtest_bits(a, state, len1, exp_bits, ctx);
            fq_nmod_mpoly_mul(b, a, t1, ctx);
            fq_nmod_mpoly_mul(t2, a, t2, ctx);

            fq_nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);

            gcd_check(g, t2, b, a, ctx, i, j, "monomial cofactors");
        }

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t1, ctx);
        fq_nmod_mpoly_clear(t2, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* one input divides the other */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_t c;
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, t1, t2;
        slong len, len1, len2;
        mp_limb_t exp_bound, exp_bound1, exp_bound2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);
        fq_nmod_init(c, ctx->fqctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t1, ctx);
        fq_nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 5);
        len1 = n_randint(state, 5);
        len2 = n_randint(state, 5);

        exp_bound = n_randint(state, 100) + 2;
        exp_bound1 = n_randint(state, 100) + 2;
        exp_bound2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(t1, state, len1, exp_bound1, ctx);
            fq_nmod_mpoly_randtest_bound(t2, state, len2, exp_bound2, ctx);
            fq_nmod_mpoly_mul(b, t1, t2, ctx);
            fq_nmod_randtest_not_zero(c, state, ctx->fqctx);
            fq_nmod_mpoly_scalar_mul_fq_nmod(a, t2, c, ctx);
            fq_nmod_randtest_not_zero(c, state, ctx->fqctx);
            fq_nmod_mpoly_scalar_mul_fq_nmod(b, b, c, ctx);

            fq_nmod_mpoly_randtest_bound(g, state, len, exp_bound, ctx);

            if ((j%2) == 0)
                fq_nmod_mpoly_swap(a, b, ctx);

            gcd_check(g, a, b, t2, ctx, i, j, "one input divides the other");
        }

        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t1, ctx);
        fq_nmod_mpoly_clear(t2, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, t;
        slong len, len1, len2;
        slong degbound;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 7, FLINT_BITS, 4);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);

        len = n_randint(state, 15) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 30/(2*ctx->minfo->nvars - 1);

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            fq_nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            fq_nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);

            fq_nmod_mpoly_mul(a, a, t, ctx);
            fq_nmod_mpoly_mul(b, b, t, ctx);

            fq_nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            gcd_check(g, a, b, t, ctx, i, j, "sparse inputs");
        }

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs with random repackings */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, t;
        mp_limb_t rlimb;
        flint_bitcnt_t newbits;
        slong len, len1, len2;
        slong degbound;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 7, FLINT_BITS, 4);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);

        len = n_randint(state, 15) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 30/(2*ctx->minfo->nvars - 1);

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            fq_nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            fq_nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
            fq_nmod_mpoly_mul(b, b, t, ctx);

            rlimb = n_randlimb(state);

            if (rlimb & UWORD(3))
            {
                newbits = a->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                fq_nmod_mpoly_repack_bits(a, a, newbits, ctx);
            }

            if (rlimb & UWORD(12))
            {
                newbits = b->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                fq_nmod_mpoly_repack_bits(b, b, newbits, ctx);
            }

            fq_nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            gcd_check(g, a, b, t, ctx, i, j, "sparse input with repacking");
        }

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs with random inflations */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, t;
        fmpz * shifts1, * shifts2, * strides;
        flint_bitcnt_t stride_bits, shift_bits;
        slong len, len1, len2;
        slong degbound;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 7, FLINT_BITS, 4);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);

        len = n_randint(state, 15) + 1;
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        degbound = 30/(2*ctx->minfo->nvars - 1);

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
            fq_nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            fq_nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            fq_nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
            fq_nmod_mpoly_mul(b, b, t, ctx);

            fq_nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts1 + k, state, shift_bits);
                fmpz_randtest_unsigned(shifts2 + k, state, shift_bits);
                fmpz_randtest_unsigned(strides + k, state, stride_bits);
            }
            fq_nmod_mpoly_inflate(a, a, shifts1, strides, ctx);
            fq_nmod_mpoly_inflate(b, b, shifts2, strides, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                if (fmpz_cmp(shifts1 + k, shifts2 + k) > 0)
                    fmpz_set(shifts1 + k, shifts2 + k);
            }

            fq_nmod_mpoly_inflate(t, t, shifts1, strides, ctx);

            gcd_check(g, a, b, t, ctx, i, j, "sparse input with inflation");
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

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* dense inputs */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, t;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 4, FLINT_BITS, 4);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);

        len1 = n_randint(state, 200) + 1;
        len2 = n_randint(state, 200);
        len3 = n_randint(state, 200);
        len4 = n_randint(state, 200);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds1[j] = 2 + n_randint(state, 12/ctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 12/ctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 12/ctx->minfo->nvars);
        }

        bits4 = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bounds(t, state, len1, degbounds1, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            fq_nmod_mpoly_randtest_bounds(a, state, len2, degbounds2, ctx);
            fq_nmod_mpoly_randtest_bounds(b, state, len3, degbounds3, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
            fq_nmod_mpoly_mul(b, b, t, ctx);

            fq_nmod_mpoly_randtest_bits(g, state, len4, bits4, ctx);

            gcd_check(g, a, b, t, ctx, i, j, "dense input");
        }

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* dense inputs with repacking */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, t;
        mp_limb_t rlimb;
        flint_bitcnt_t newbits;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 4, FLINT_BITS, 4);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);

        len1 = n_randint(state, 200) + 1;
        len2 = n_randint(state, 200);
        len3 = n_randint(state, 200);
        len4 = n_randint(state, 200);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds1[j] = 2 + n_randint(state, 16/ctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 16/ctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 16/ctx->minfo->nvars);
        }

        bits4 = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bounds(t, state, len1, degbounds1, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            fq_nmod_mpoly_randtest_bounds(a, state, len2, degbounds2, ctx);
            fq_nmod_mpoly_randtest_bounds(b, state, len3, degbounds3, ctx);

            fq_nmod_mpoly_mul(a, a, t, ctx);
            fq_nmod_mpoly_mul(b, b, t, ctx);

            rlimb = n_randlimb(state);

            if (rlimb & UWORD(3))
            {
                newbits = a->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                fq_nmod_mpoly_repack_bits(a, a, newbits, ctx);
            }

            if (rlimb & UWORD(12))
            {
                newbits = b->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                fq_nmod_mpoly_repack_bits(b, b, newbits, ctx);
            }

            fq_nmod_mpoly_randtest_bits(g, state, len4, bits4, ctx);

            gcd_check(g, a, b, t, ctx, i, j, "dense input with repacking");
        }

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* dense inputs with random inflations */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t a, b, g, t;
        fmpz * shifts1, * shifts2, * strides;
        flint_bitcnt_t stride_bits, shift_bits;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 4, FLINT_BITS, 4);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(a, ctx);
        fq_nmod_mpoly_init(b, ctx);
        fq_nmod_mpoly_init(t, ctx);

        len1 = n_randint(state, 200) + 1;
        len2 = n_randint(state, 200);
        len3 = n_randint(state, 200);
        len4 = n_randint(state, 200);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds1[j] = 2 + n_randint(state, 15/ctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 15/ctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 15/ctx->minfo->nvars);
        }

        bits4 = n_randint(state, 200);

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
            fq_nmod_mpoly_randtest_bounds(t, state, len1, degbounds1, ctx);
            if (fq_nmod_mpoly_is_zero(t, ctx))
                fq_nmod_mpoly_one(t, ctx);
            fq_nmod_mpoly_randtest_bounds(a, state, len2, degbounds2, ctx);
            fq_nmod_mpoly_randtest_bounds(b, state, len3, degbounds3, ctx);
            fq_nmod_mpoly_mul(a, a, t, ctx);
            fq_nmod_mpoly_mul(b, b, t, ctx);

            fq_nmod_mpoly_randtest_bits(g, state, len4, bits4, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts1 + k, state, shift_bits);
                fmpz_randtest_unsigned(shifts2 + k, state, shift_bits);
                fmpz_randtest_unsigned(strides + k, state, stride_bits);
            }
            fq_nmod_mpoly_inflate(a, a, shifts1, strides, ctx);
            fq_nmod_mpoly_inflate(b, b, shifts2, strides, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                if (fmpz_cmp(shifts1 + k, shifts2 + k) > 0)
                    fmpz_set(shifts1 + k, shifts2 + k);
            }
            fq_nmod_mpoly_inflate(t, t, shifts1, strides, ctx);

            gcd_check(g, a, b, t, ctx, i, j, "dense input with inflation");
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

        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(a, ctx);
        fq_nmod_mpoly_clear(b, ctx);
        fq_nmod_mpoly_clear(t, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef gcd_check
