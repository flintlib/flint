/*
    Copyright (C) 2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_divides_monagan_pearce, state)
{
    int i, j, result, ok1, ok2;

    /*
        A bad case is hit when testing with multiplier 50. The following
        example illustrates this behaviour if the ordering is changed to
        ORD_DEGLEX
    */

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, q, r;

        fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(q, ctx);
        fmpz_mpoly_init(r, ctx);

        fmpz_mpoly_set_str_pretty(f, "-x1^1918*x2^1075-x1^1891*x2^2001",NULL, ctx);
        fmpz_mpoly_set_str_pretty(g, "x1^22*x2^3-x1^19*x2^21-x1^16*x2^10-2*x1^14*x2^17-x1^14*x2^11-x1*x2^15-2*x2^17", NULL, ctx);

        ok1 = fmpz_mpoly_divides_monagan_pearce(q, f, g, ctx);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(q, ctx);
        fmpz_mpoly_clear(r, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10) + 1;

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            do {
                fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits + 1, exp_bits2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            ok1 = fmpz_mpoly_divides_monagan_pearce(k, h, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);
            result = (ok1 && fmpz_mpoly_equal(f, k, ctx));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check random polys don't divide */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits;
        fmpz * shifts, * strides;
        slong n;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound =  n_randint(state, 1000/n/n) + 1;
        exp_bound1 = n_randint(state, 1000/n/n) + 1;
        exp_bound2 = n_randint(state, 1000/n/n) + 1;

        coeff_bits = n_randint(state, 200);

        shifts = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        strides = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            fmpz_init(shifts + j);
            fmpz_init(strides + j);
            fmpz_randtest_unsigned(shifts + j, state, 100);
            fmpz_randtest_unsigned(strides + j, state, 100);
            fmpz_add_ui(strides + j, strides + j, 1);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_inflate(f, f, shifts, strides, ctx);
            fmpz_mpoly_inflate(g, g, shifts, strides, ctx);

            ok1 = fmpz_mpoly_divides_monagan_pearce(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);

            if (ok1 == 0)
                continue;

            fmpz_mpoly_mul_johnson(k, h, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);
            result = fmpz_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check random polys don't divide\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            fmpz_clear(shifts + j);
            fmpz_clear(strides + j);
        }
        flint_free(shifts);
        flint_free(strides);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument, exact division */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50) + 1;

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            do {
                fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits + 1, exp_bits2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            ok1 = fmpz_mpoly_divides_monagan_pearce(k, h, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);
            ok2 = fmpz_mpoly_divides_monagan_pearce(h, h, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            result = (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(h, k, ctx));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first argument, exact division\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing, first argument, random polys */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;
       slong n;

       fmpz_mpoly_ctx_init_rand(ctx, state, 10);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       len = n_randint(state, 20);
       len1 = n_randint(state, 20);
       len2 = n_randint(state, 20) + 1;

       n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
       exp_bound =  n_randint(state, 1000/n/n) + 1;
       exp_bound1 = n_randint(state, 1000/n/n) + 1;
       exp_bound2 = n_randint(state, 1000/n/n) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

          ok1 = fmpz_mpoly_divides_monagan_pearce(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          ok2 = fmpz_mpoly_divides_monagan_pearce(f, f, g, ctx);
          fmpz_mpoly_assert_canonical(f, ctx);

          result = ((ok1 == ok2) &&  (ok1 == 0 || fmpz_mpoly_equal(f, h, ctx)));

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing, first argument, random polys\ni = %wd, j = %wd\n", i ,j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument, exact division */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50) + 1;

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            do {
                fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits + 1, exp_bits2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            ok1 = fmpz_mpoly_divides_monagan_pearce(k, h, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);
            ok2 = fmpz_mpoly_divides_monagan_pearce(g, h, g, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);
            result = (ok1 == 1 && ok2 == 1 && fmpz_mpoly_equal(g, k, ctx));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second argument, exact division\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing, second argument, random polys */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits;
       slong n;

       fmpz_mpoly_ctx_init_rand(ctx, state, 10);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100) + 1;

       n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
       exp_bound =  n_randint(state, 1000/n/n) + 1;
       exp_bound1 = n_randint(state, 1000/n/n) + 1;
       exp_bound2 = n_randint(state, 1000/n/n) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          do {
             fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

          ok1 = fmpz_mpoly_divides_monagan_pearce(h, f, g, ctx);
          fmpz_mpoly_assert_canonical(h, ctx);
          ok2 = fmpz_mpoly_divides_monagan_pearce(g, f, g, ctx);
          fmpz_mpoly_assert_canonical(g, ctx);

          result = ((ok1 == ok2) &&  (ok1 == 0 || fmpz_mpoly_equal(g, h, ctx)));

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing, second argument, random polys\ni = %wd, j = %wd\n", i ,j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
