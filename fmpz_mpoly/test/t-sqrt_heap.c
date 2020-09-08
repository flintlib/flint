/*
    Copyright (C) 2017, 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mpoly.h"

int
main(void)
{
    slong i, j, tmul = 10;
    FLINT_TEST_INIT(state);

    flint_printf("sqrt_heap....");
    fflush(stdout);

    {
        fmpz_mpoly_t f, g, p, q;
        fmpz_mpoly_ctx_t ctx;
        const char * vars[] = {"x", "y", "z", "t", "u"};
        int sqr;

        fmpz_mpoly_ctx_init(ctx, 5, ORD_LEX);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(p, ctx);
        fmpz_mpoly_init(q, ctx);

        fmpz_mpoly_set_str_pretty(f, "(1+x+y+2*z^2+3*t^3+5*u^5)^6", vars, ctx);

        fmpz_mpoly_mul(p, f, f, ctx);
        fmpz_mpoly_assert_canonical(p, ctx);

        sqr = fmpz_mpoly_sqrt_heap(q, p, ctx, 1);
        fmpz_mpoly_assert_canonical(q, ctx);

	     if (!sqr)
	     {
	         printf("FAIL\n");
            flint_printf("Check example1: sqr\n");
	         abort();
	     }

        fmpz_mpoly_mul(g, q, q, ctx);

        if (!fmpz_mpoly_equal(p, g, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check example1\n");
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(p, ctx);
        fmpz_mpoly_clear(q, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check sqrt(f^2)^2 = f^2 */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1;
        flint_bitcnt_t coeff_bits, coeff_bits1;
        int sqr;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);

        exp_bits =  n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 200);
        coeff_bits1 = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits1, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_mul(g, f, f, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

   	      sqr = fmpz_mpoly_sqrt_heap(h, g, ctx, 1);
            fmpz_mpoly_assert_canonical(h, ctx);

            if (!sqr)
            {
                printf("FAIL\n");
                flint_printf("Check sqrt(f^2) returns 1\n");
                flint_abort();
            }

            fmpz_mpoly_mul(k, h, h, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);

            if (!fmpz_mpoly_equal(g, k, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check sqrt(f^2)^2 = f^2\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            sqr = fmpz_mpoly_sqrt_heap(h, g, ctx, 0);
            fmpz_mpoly_assert_canonical(h, ctx);

            if (!sqr)
            {
                printf("FAIL\n");
                flint_printf("Check sqrt(f^2) returns 1: nocheck\n");
                flint_abort();
            }

            fmpz_mpoly_mul(k, h, h, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);

            if (!fmpz_mpoly_equal(g, k, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check sqrt(f^2)^2 = f^2\ni = %wd, j = %wd: nocheck\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
    }

    /* Check sqrt(f^2*(x^2+x)) returns 0 */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k, x;
        slong len, len1, nvars;
        flint_bitcnt_t exp_bits, exp_bits1;
        flint_bitcnt_t coeff_bits, coeff_bits1;
        int sqr;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);
        fmpz_mpoly_init(x, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100) + 1;

        exp_bits =  n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 200);
        coeff_bits1 = n_randint(state, 200) + 1;

        nvars = fmpz_mpoly_ctx_nvars(ctx);

        for (j = 0; j < 4; j++)
        {
            do {
               fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits1, exp_bits1, ctx);
             } while (fmpz_mpoly_is_zero(f, ctx));
            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_mul(g, f, f, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

            fmpz_mpoly_gen(x, n_randint(state, nvars), ctx);
            fmpz_mpoly_mul(k, x, x, ctx);
            fmpz_mpoly_add(k, k, x, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);

            fmpz_mpoly_mul(g, g, k, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

   	      sqr = fmpz_mpoly_sqrt_heap(h, g, ctx, 1);

            fmpz_mpoly_assert_canonical(h, ctx);

            if (sqr)
            {
                printf("FAIL\n");
                flint_printf("Check nonsquare returns 0\n");
                flint_abort();
            }

            if (!fmpz_mpoly_is_zero(h, ctx))
            {
               printf("FAIL\n");
               flint_printf("Nonsquare returns 0 sqrt\n");
               flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_clear(x, ctx);
    }

    /* Check aliasing of square root with input */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1;
        flint_bitcnt_t coeff_bits, coeff_bits1;
        int sqr1, sqr2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);

        exp_bits =  n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 200);
        coeff_bits1 = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits1, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_mul(g, f, f, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);
            fmpz_mpoly_set(k, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);

            sqr1 = fmpz_mpoly_sqrt_heap(h, g, ctx, 1);
            fmpz_mpoly_assert_canonical(h, ctx);

            sqr2 = fmpz_mpoly_sqrt_heap(g, g, ctx, 1);
            fmpz_mpoly_assert_canonical(g, ctx);

            if (sqr1 != sqr2 || !fmpz_mpoly_equal(g, h, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check aliasing\n");
                flint_abort();
            }

            sqr2 = fmpz_mpoly_sqrt_heap(k, k, ctx, 0);
            fmpz_mpoly_assert_canonical(k, ctx);

            if (sqr1 != sqr2 || !fmpz_mpoly_equal(k, h, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check aliasing: nocheck\n");
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

