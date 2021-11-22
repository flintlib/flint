/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mod_mpoly.h"

int
main(void)
{
    slong i, j, w;
    FLINT_TEST_INIT(state);

    flint_printf("divrem_ideal_monagan_pearce....");
    fflush(stdout);

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h, k, r;
        slong len, len1, len2;
        slong exp_bits, exp_bits1, exp_bits2;
        fmpz_mod_mpoly_struct * qarr[1], * darr[1];

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 20, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_mod_mpoly_init(k, ctx);
        fmpz_mod_mpoly_init(r, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100) + 1;

        exp_bits = n_randint(state, 100) + 1;
        exp_bits1 = n_randint(state, 100) + 1;
        exp_bits2 = n_randint(state, 100) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len2, exp_bits2 + 1, ctx);
            if (fmpz_mod_mpoly_is_zero(g, ctx))
                fmpz_mod_mpoly_one(g, ctx);
            fmpz_mod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_randtest_bits(r, state, len, exp_bits, ctx);

            fmpz_mod_mpoly_mul(h, f, g, ctx);

            qarr[0] = k;
            darr[0] = g;

            fmpz_mod_mpoly_divrem_ideal_monagan_pearce(qarr, r, h, darr, 1, ctx);
            fmpz_mod_mpoly_assert_canonical(qarr[0], ctx);
            fmpz_mod_mpoly_assert_canonical(r, ctx);

            if (!fmpz_mod_mpoly_equal(f, k, ctx))
            {
                flint_printf("FAIL: Check f*g/g = f\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_clear(k, ctx);
        fmpz_mod_mpoly_clear(r, ctx);

        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check f = g1*q1 + ... + gn*qn + r for random polys */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, r, k1, k2;
        fmpz_mod_mpoly_struct * g, * q;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
        fmpz_mod_mpoly_struct * qarr[5], * darr[5];
        slong n;

        num = n_randint(state, 5) + 1;

        g = FLINT_ARRAY_ALLOC(num, fmpz_mod_mpoly_struct);
        q = FLINT_ARRAY_ALLOC(num, fmpz_mod_mpoly_struct);

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 10, 200);

        for (w = 0; w < num; w++)
        {
            fmpz_mod_mpoly_init(g + w, ctx);
            darr[w] = g + w;

            fmpz_mod_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }  

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(k1, ctx);
        fmpz_mod_mpoly_init(k2, ctx);
        fmpz_mod_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound = n_randint(state, 10/n + 1) + 2;
        exp_bound1 = n_randint(state, 25/n + 1) + 2;
        exp_bound2 = n_randint(state, 20/n + 1) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_randtest_bound(darr[w], state, len2, exp_bound2 + 1, ctx);
                if (fmpz_mod_mpoly_is_zero(darr[w], ctx))
                    fmpz_mod_mpoly_one(darr[w], ctx);
                fmpz_mod_mpoly_randtest_bound(qarr[w], state, len, exp_bound, ctx);
            }
            fmpz_mod_mpoly_randtest_bound(k1, state, len, exp_bound, ctx);
            fmpz_mod_mpoly_randtest_bound(k2, state, len, exp_bound, ctx);

            fmpz_mod_mpoly_divrem_ideal_monagan_pearce(qarr, r, f, darr, num, ctx);
            fmpz_mod_mpoly_assert_canonical(r, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_assert_canonical(qarr[w], ctx);
                fmpz_mod_mpoly_remainder_strongtest(r, darr[w], ctx);
            }

            fmpz_mod_mpoly_zero(k2, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
                fmpz_mod_mpoly_add(k2, k2, k1, ctx);
            }
            fmpz_mod_mpoly_add(k2, k2, r, ctx);

            if (!fmpz_mod_mpoly_equal(f, k2, ctx))
            {
                flint_printf("FAIL: Check f = g1*q1 + ... + gn*qn + r\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (w = 0; w < num; w++)
            fmpz_mod_mpoly_clear(qarr[w], ctx);
        for (w = 0; w < num; w++)
            fmpz_mod_mpoly_clear(darr[w], ctx);
        fmpz_mod_mpoly_clear(f, ctx);  
        fmpz_mod_mpoly_clear(k1, ctx);  
        fmpz_mod_mpoly_clear(k2, ctx);  
        fmpz_mod_mpoly_clear(r, ctx);

        flint_free(g);
        flint_free(q);  
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, r, k1, k2;
        fmpz_mod_mpoly_struct * g, * q;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
        fmpz_mod_mpoly_struct * qarr[5], * darr[5];
        slong n;

        num = n_randint(state, 5) + 1;

        g = FLINT_ARRAY_ALLOC(num, fmpz_mod_mpoly_struct);
        q = FLINT_ARRAY_ALLOC(num, fmpz_mod_mpoly_struct);

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 10, 200);

        for (w = 0; w < num; w++)
        {
            fmpz_mod_mpoly_init(g + w, ctx);
            darr[w] = g + w;

            fmpz_mod_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }  

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(k1, ctx);
        fmpz_mod_mpoly_init(k2, ctx);
        fmpz_mod_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound = n_randint(state, 10/n + 1) + 2;
        exp_bound1 = n_randint(state, 25/n + 1) + 2;
        exp_bound2 = n_randint(state, 20/n + 1) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_randtest_bound(darr[w], state, len2, exp_bound2 + 1, ctx);
                if (fmpz_mod_mpoly_is_zero(darr[w], ctx))
                    fmpz_mod_mpoly_one(darr[w], ctx);
                fmpz_mod_mpoly_randtest_bound(qarr[w], state, len, exp_bound, ctx);
            }
            fmpz_mod_mpoly_randtest_bound(k1, state, len, exp_bound, ctx);
            fmpz_mod_mpoly_randtest_bound(k2, state, len, exp_bound, ctx);

            fmpz_mod_mpoly_set(r, f, ctx);

            fmpz_mod_mpoly_divrem_ideal_monagan_pearce(qarr, f, f, darr, num, ctx);
            fmpz_mod_mpoly_assert_canonical(f, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_assert_canonical(qarr[w], ctx);
                fmpz_mod_mpoly_remainder_strongtest(f, darr[w], ctx);
            }

            fmpz_mod_mpoly_zero(k2, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mod_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
                fmpz_mod_mpoly_add(k2, k2, k1, ctx);
            }
            fmpz_mod_mpoly_add(k2, k2, f, ctx);

            if (!fmpz_mod_mpoly_equal(r, k2, ctx))
            {
                flint_printf("FAIL: Check aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (w = 0; w < num; w++)
            fmpz_mod_mpoly_clear(qarr[w], ctx);
        for (w = 0; w < num; w++)
            fmpz_mod_mpoly_clear(darr[w], ctx);
        fmpz_mod_mpoly_clear(f, ctx);  
        fmpz_mod_mpoly_clear(k1, ctx);  
        fmpz_mod_mpoly_clear(k2, ctx);  
        fmpz_mod_mpoly_clear(r, ctx);

        flint_free(g);
        flint_free(q);  
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    flint_printf("PASS\n");
    FLINT_TEST_CLEANUP(state);

    return 0;
}

