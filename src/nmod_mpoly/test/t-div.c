/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_div, state)
{
    slong i, j, tmul = 10;

    {
        nmod_mpoly_t f, g, p, q;
        nmod_mpoly_ctx_t ctx;
        const char * vars[] = {"x", "y", "z", "t", "u"};

        nmod_mpoly_ctx_init(ctx, 5, ORD_LEX, 1000003);
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(p, ctx);
        nmod_mpoly_init(q, ctx);

        nmod_mpoly_set_str_pretty(f, "(1+x+y+2*z^2+3*t^3+5*u^5)^6", vars, ctx);
        nmod_mpoly_set_str_pretty(g, "(1+u+t+2*z^2+3*y^3+5*x^5)^6", vars, ctx);

        nmod_mpoly_mul(p, f, g, ctx);
        nmod_mpoly_assert_canonical(p, ctx);

        nmod_mpoly_div(q, p, f, ctx);
        nmod_mpoly_assert_canonical(q, ctx);

        if (!nmod_mpoly_equal(q, g, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check example\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_clear(q, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f */
    for (i = 0; i < 10 * tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k, l;
        mp_limb_t modulus;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);
        nmod_mpoly_init(l, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50) + 1;

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bits2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            nmod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);
            nmod_mpoly_randtest_bits(l, state, len, exp_bits, ctx);

            nmod_mpoly_mul(h, f, g, ctx);

            nmod_mpoly_div(k, h, g, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            if (!nmod_mpoly_equal(k, f, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_set(l, h, ctx);
            nmod_mpoly_div(l, l, g, ctx);
            nmod_mpoly_assert_canonical(l, ctx);
            if (!nmod_mpoly_equal(l, f, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f aliasing dividend\n"
                                                   "i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_div(g, h, g, ctx);
            nmod_mpoly_assert_canonical(g, ctx);
            if (!nmod_mpoly_equal(g, f, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f aliasing divisor\n"
                                                   "i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_clear(l, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check div matches divrem for random polys */
    for (i = 0; i < 5 * tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, q, r, k;
        mp_limb_t modulus;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        fmpz * shifts, * strides;
        slong n, nvars;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);
        nvars = ctx->minfo->nvars;

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(q, ctx);
        nmod_mpoly_init(r, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), nvars);
        max_bound = 1 + 400/n/n;
        exp_bound = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        shifts = (fmpz *) flint_malloc(nvars*sizeof(fmpz));
        strides = (fmpz *) flint_malloc(nvars*sizeof(fmpz));
        for (j = 0; j < nvars; j++)
        {
            exp_bound[j] = UWORD(1) << (FLINT_BITS - 1);
            exp_bound1[j] = n_randint(state, max_bound) + 1;
            exp_bound2[j] = n_randint(state, max_bound) + 1;
            fmpz_init(shifts + j);
            fmpz_init(strides + j);
            fmpz_randtest_unsigned(shifts + j, state, 100);
            fmpz_randtest_unsigned(strides + j, state, 100);
            fmpz_add_ui(strides + j, strides + j, 1);
        }

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bounds(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bounds(g, state, len2, exp_bound2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bounds(q, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bounds(r, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bounds(k, state, len, exp_bound, ctx);

            nmod_mpoly_inflate(f, f, shifts, strides, ctx);
            nmod_mpoly_inflate(g, g, shifts, strides, ctx);

            nmod_mpoly_divrem(q, r, f, g, ctx);
            nmod_mpoly_assert_canonical(q, ctx);
            nmod_mpoly_assert_canonical(r, ctx);
            nmod_mpoly_remainder_strongtest(r, g, ctx);

            nmod_mpoly_mul(k, q, g, ctx);
            nmod_mpoly_add(k, k, r, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            if (!nmod_mpoly_equal(f, k, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check f = g*q + r for random polys\n"
                                                   "i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_div(k, f, g, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            if (!nmod_mpoly_equal(k, q, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check div matches divrem\n"
                                                   "i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_set(k, f, ctx);
            nmod_mpoly_div(k, k, g, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            if (!nmod_mpoly_equal(k, q, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check div matches divrem aliasing dividend\n"
                                                   "i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_div(g, f, g, ctx);
            nmod_mpoly_assert_canonical(g, ctx);
            if (!nmod_mpoly_equal(g, q, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check div matches divrem aliasing divisor\n"
                                                   "i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (j = 0; j < nvars; j++)
        {
            fmpz_clear(shifts + j);
            fmpz_clear(strides + j);
        }
        flint_free(shifts);
        flint_free(strides);

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(q, ctx);
        nmod_mpoly_clear(r, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
