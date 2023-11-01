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
#define gcd_check gcd_check_gcd
void gcd_check(
    nmod_mpoly_t g,
    nmod_mpoly_t a,
    nmod_mpoly_t b,
    nmod_mpoly_ctx_t ctx,
    slong i,
    slong j,
    const char * name)
{
    nmod_mpoly_t ca, cb, cg;

    nmod_mpoly_init(ca, ctx);
    nmod_mpoly_init(cb, ctx);
    nmod_mpoly_init(cg, ctx);

    if (!nmod_mpoly_gcd(g, a, b, ctx))
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

    if (((i + j) % 7) == 0)
    {
        nmod_mpoly_set(cg, b, ctx);
        nmod_mpoly_gcd(cg, cg, a, ctx);
        if (!nmod_mpoly_equal(cg, g, ctx))
        {
            flint_printf("FAIL: check aliasing 1\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    if (((i + j) % 9) == 0)
    {
        nmod_mpoly_set(cg, b, ctx);
        nmod_mpoly_gcd(cg, a, cg, ctx);
        if (!nmod_mpoly_equal(cg, g, ctx))
        {
            flint_printf("FAIL: check aliasing 2\n");
            flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
            fflush(stdout);
            flint_abort();
        }
    }

    if (!nmod_mpoly_divides(ca, a, g, ctx) ||
        !nmod_mpoly_divides(cb, b, g, ctx))
    {
        flint_printf("FAIL: check divisibility\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    if (!nmod_mpoly_gcd(cg, ca, cb, ctx))
    {
        flint_printf("FAIL: check cofactor gcd can be computed\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

    nmod_mpoly_assert_canonical(cg, ctx);

    if (!nmod_mpoly_is_one(cg, ctx))
    {
        flint_printf("FAIL: check cofactor gcd is one\n");
        flint_printf("i = %wd, j = %wd, %s\n", i, j, name);
        fflush(stdout);
        flint_abort();
    }

cleanup:

    nmod_mpoly_clear(ca, ctx);
    nmod_mpoly_clear(cb, ctx);
    nmod_mpoly_clear(cg, ctx);
}

TEST_FUNCTION_START(nmod_mpoly_gcd, state)
{
    const slong max_threads = 5;
    slong i, j, k, tmul = 5;

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t g, a, b;
        const char * vars[] = {"t" ,"x", "y", "z"};

        nmod_mpoly_ctx_init(ctx, 4, ORD_LEX, 536870909);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(g, ctx);

        nmod_mpoly_set_str_pretty(g, "x + y + z + t", vars, ctx);
        nmod_mpoly_set_str_pretty(a, "x^2 + y^2 + z^2 + t^2", vars, ctx);
        nmod_mpoly_set_str_pretty(b, "x^3 + y^3 + z^3 + t^3", vars, ctx);
        nmod_mpoly_mul(a, a, g, ctx);
        nmod_mpoly_mul(b, b, g, ctx);

        gcd_check(g, a, b, ctx, 0, 0, "sparse example");

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t g, a, b;
        const char * vars[] = {"t" ,"x", "y", "z"};

        nmod_mpoly_ctx_init(ctx, 4, ORD_LEX, 536870909);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(g, ctx);

        nmod_mpoly_set_str_pretty(g, "39 - t*x + 39*x^100 - t*x^101 + 39*x^3*y - t*x^4*y - 7*x^2*y^3*z^11 - 7*x^102*y^3*z^11 - 7*x^5*y^4*z^11 + 78*t^15*x^78*y^3*z^13 - 2*t^16*x^79*y^3*z^13 + x^1000*y^3*z^20 + x^1100*y^3*z^20 + x^1003*y^4*z^20 - 14*t^15*x^80*y^6*z^24 + 2*t^15*x^1078*y^6*z^33", vars, ctx);
        nmod_mpoly_set_str_pretty(a, "39 - t*x - 7*x^2*y^3*z^11 + x^1000*y^3*z^20", vars, ctx);
        nmod_mpoly_set_str_pretty(b, "1 + x^100 + x^3*y + 2*t^15*x^78*y^3*z^13", vars, ctx);
        nmod_mpoly_mul(a, a, g, ctx);
        nmod_mpoly_mul(b, b, g, ctx);

        gcd_check(g, a, b, ctx, 0, 0, "sparse example");

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t g, a, b;
        const char * vars[] = {"x", "y"};

        nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, 536870909);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(g, ctx);

        nmod_mpoly_set_str_pretty(g, "(1+x+y)^100", vars, ctx);
        nmod_mpoly_set_str_pretty(a, "(2+x+y)^100", vars, ctx);
        nmod_mpoly_set_str_pretty(b, "(3+x+y)^100", vars, ctx);
        nmod_mpoly_mul(a, a, g, ctx);
        nmod_mpoly_mul(b, b, g, ctx);

        gcd_check(g, a, b, ctx, 0, 0, "big example");

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    for (i = 3; i <= 8; i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t g, a, b, t;

        nmod_mpoly_ctx_init(ctx, i, ORD_DEGREVLEX, 43051);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(t, ctx);

        nmod_mpoly_one(g, ctx);
        nmod_mpoly_one(a, ctx);
        nmod_mpoly_one(b, ctx);
        for (j = 0; j < i; j++)
        {
            nmod_mpoly_gen(t, j, ctx);
            nmod_mpoly_add_ui(t, t, 1, ctx);
            nmod_mpoly_mul(g, g, t, ctx);
            nmod_mpoly_gen(t, j, ctx);
            nmod_mpoly_sub_ui(t, t, 2, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_gen(t, j, ctx);
            nmod_mpoly_add_ui(t, t, 2, ctx);
            nmod_mpoly_mul(b, b, t, ctx);
        }
        nmod_mpoly_sub_ui(g, g, 2, ctx);
        nmod_mpoly_add_ui(a, a, 2, ctx);
        nmod_mpoly_sub_ui(b, b, 2, ctx);

        nmod_mpoly_mul(a, a, g, ctx);
        nmod_mpoly_mul(b, b, g, ctx);

        gcd_check(g, a, b, ctx, i, 0, "dense examples");

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t g, a, b;
        const char * vars[] = {"x1", "x2", "x3", "x4"};

        nmod_mpoly_ctx_init(ctx, 2, ORD_DEGREVLEX, 2);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(g, ctx);

        nmod_mpoly_set_str_pretty(a, "x1*x2^5+x1*x2^3+x2^4+x2^2", vars, ctx);
        nmod_mpoly_set_str_pretty(b, "x1*x2^7+x1*x2^6+x1*x2^5+x2^6+x2^5+x2^4+x1*x2^2+x2", vars, ctx);
        nmod_mpoly_set_str_pretty(g, "1", vars, ctx);
        nmod_mpoly_mul(a, a, g,ctx);
        nmod_mpoly_mul(b, b, g,ctx);

        flint_set_num_threads(n_randint(state, max_threads) + 1);
        gcd_check(g, a, b, ctx, 0, 0, "example");

        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* The gcd should always work when one input is a monomial */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
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
            nmod_mpoly_randtest_bits(t, state, 1, exp_bits, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            nmod_mpoly_randtest_bits(a, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(b, state, len2, exp_bits2, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_mul(b, b, t, ctx);

            nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "monomial");
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* The gcd should always work when both cofactors are monomials */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t1, t2;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t1, ctx);
        nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 25);
        len1 = n_randint(state, 25);

        exp_bits = n_randint(state, 70) + 2;
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(t1, state, 1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(t2, state, 1, exp_bits2, ctx);
            if (t1->length != 1 || t2->length != 1)
            {
                flint_printf("FAIL:\ncheck random monomial generation\n");
                fflush(stdout);
                flint_abort();
            }
            nmod_mpoly_randtest_bits(a, state, len1, exp_bits, ctx);
            nmod_mpoly_mul(b, a, t1, ctx);
            nmod_mpoly_mul(a, a, t2, ctx);

            nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "monomial cofactors");
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t1, ctx);
        nmod_mpoly_clear(t2, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* one input divides the other */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        mp_limb_t c;
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t1, t2;
        slong len, len1, len2;
        mp_limb_t exp_bound, exp_bound1, exp_bound2;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t1, ctx);
        nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        exp_bound = n_randint(state, 100) + 2;
        exp_bound1 = n_randint(state, 100) + 2;
        exp_bound2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(t1, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(t2, state, len2, exp_bound2, ctx);
            nmod_mpoly_mul(b, t1, t2, ctx);
            c = n_randint(state, ctx->mod.n);
            nmod_mpoly_scalar_mul_ui(a, t2, c, ctx);
            c = n_randint(state, ctx->mod.n);
            nmod_mpoly_scalar_mul_ui(b, b, c, ctx);

            nmod_mpoly_randtest_bound(g, state, len, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            if ((j%2) == 0)
                nmod_mpoly_swap(a, b, ctx);

            gcd_check(g, a, b, ctx, i, j, "one input divides the other");
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t1, ctx);
        nmod_mpoly_clear(t2, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t;
        slong len, len1, len2;
        slong degbound;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 7, modulus);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);

        degbound = 40/(2*ctx->minfo->nvars - 1);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);

            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_mul(b, b, t, ctx);

            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "sparse inputs");
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs with random repackings */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t;
        mp_limb_t rlimb;
        flint_bitcnt_t newbits;
        slong len, len1, len2;
        slong degbound;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 7, modulus);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);

        degbound = 40/(2*ctx->minfo->nvars - 1);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_mul(b, b, t, ctx);

            rlimb = n_randlimb(state);

            if (rlimb & UWORD(3))
            {
                newbits = a->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                nmod_mpoly_repack_bits(a, a, newbits, ctx);
            }

            if (rlimb & UWORD(12))
            {
                newbits = b->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                nmod_mpoly_repack_bits(b, b, newbits, ctx);
            }

            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "sparse input with repacking");
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* sparse inputs with random inflations */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t;
        fmpz * shifts1, * shifts2, * strides;
        flint_bitcnt_t stride_bits, shift_bits;
        slong len, len1, len2;
        slong degbound;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 7, modulus);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);

        len = n_randint(state, 20) + 1;
        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30);

        degbound = 40/(2*ctx->minfo->nvars - 1);

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
            nmod_mpoly_randtest_bound(t, state, len, degbound, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            nmod_mpoly_randtest_bound(a, state, len1, degbound, ctx);
            nmod_mpoly_randtest_bound(b, state, len2, degbound, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_mul(b, b, t, ctx);

            nmod_mpoly_randtest_bits(g, state, len, FLINT_BITS, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts1 + k, state, shift_bits);
                fmpz_randtest_unsigned(shifts2 + k, state, shift_bits);
                fmpz_randtest_unsigned(strides + k, state, stride_bits);
            }
            nmod_mpoly_inflate(a, a, shifts1, strides, ctx);
            nmod_mpoly_inflate(b, b, shifts2, strides, ctx);

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

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* dense inputs */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 4, modulus);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);

        len1 = n_randint(state, 300) + 1;
        len2 = n_randint(state, 300);
        len3 = n_randint(state, 300);
        len4 = n_randint(state, 300);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds1[j] = 2 + n_randint(state, 16/ctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 16/ctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 16/ctx->minfo->nvars);
        }

        bits4 = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bounds(t, state, len1, degbounds1, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            nmod_mpoly_randtest_bounds(a, state, len2, degbounds2, ctx);
            nmod_mpoly_randtest_bounds(b, state, len3, degbounds3, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_mul(b, b, t, ctx);

            nmod_mpoly_randtest_bits(g, state, len4, bits4, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "dense input");
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* dense inputs with repacking */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t;
        mp_limb_t rlimb;
        flint_bitcnt_t newbits;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 4, modulus);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);

        len1 = n_randint(state, 300) + 1;
        len2 = n_randint(state, 300);
        len3 = n_randint(state, 300);
        len4 = n_randint(state, 300);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            degbounds1[j] = 2 + n_randint(state, 16/ctx->minfo->nvars);
            degbounds2[j] = 1 + n_randint(state, 16/ctx->minfo->nvars);
            degbounds3[j] = 1 + n_randint(state, 16/ctx->minfo->nvars);
        }

        bits4 = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bounds(t, state, len1, degbounds1, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            nmod_mpoly_randtest_bounds(a, state, len2, degbounds2, ctx);
            nmod_mpoly_randtest_bounds(b, state, len3, degbounds3, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_mul(b, b, t, ctx);

            rlimb = n_randlimb(state);

            if (rlimb & UWORD(3))
            {
                newbits = a->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                nmod_mpoly_repack_bits(a, a, newbits, ctx);
            }

            if (rlimb & UWORD(12))
            {
                newbits = b->bits + n_randint(state, 2*FLINT_BITS);
                newbits = mpoly_fix_bits(newbits, ctx->minfo);
                nmod_mpoly_repack_bits(b, b, newbits, ctx);
            }

            nmod_mpoly_randtest_bits(g, state, len4, bits4, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);
            gcd_check(g, a, b, ctx, i, j, "dense input with repacking");
        }

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* dense inputs with random inflations */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t a, b, g, t;
        fmpz * shifts1, * shifts2, * strides;
        flint_bitcnt_t stride_bits, shift_bits;
        slong len1, len2, len3, len4;
        ulong degbounds1[4];
        ulong degbounds2[4];
        ulong degbounds3[4];
        flint_bitcnt_t bits4;
        mp_limb_t modulus;

        modulus = n_randint(state, (i % 10 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 4, modulus);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(a, ctx);
        nmod_mpoly_init(b, ctx);
        nmod_mpoly_init(t, ctx);

        len1 = n_randint(state, 300) + 1;
        len2 = n_randint(state, 300);
        len3 = n_randint(state, 300);
        len4 = n_randint(state, 300);

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
            nmod_mpoly_randtest_bounds(t, state, len1, degbounds1, ctx);
            if (nmod_mpoly_is_zero(t, ctx))
                nmod_mpoly_one(t, ctx);
            nmod_mpoly_randtest_bounds(a, state, len2, degbounds2, ctx);
            nmod_mpoly_randtest_bounds(b, state, len3, degbounds3, ctx);
            nmod_mpoly_mul(a, a, t, ctx);
            nmod_mpoly_mul(b, b, t, ctx);

            nmod_mpoly_randtest_bits(g, state, len4, bits4, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts1 + k, state, shift_bits);
                fmpz_randtest_unsigned(shifts2 + k, state, shift_bits);
                fmpz_randtest_unsigned(strides + k, state, stride_bits);
            }
            nmod_mpoly_inflate(a, a, shifts1, strides, ctx);
            nmod_mpoly_inflate(b, b, shifts2, strides, ctx);

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

        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(a, ctx);
        nmod_mpoly_clear(b, ctx);
        nmod_mpoly_clear(t, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#undef gcd_check
