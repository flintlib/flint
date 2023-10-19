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

TEST_FUNCTION_START(fmpz_mpoly_divrem_ideal_monagan_pearce, state)
{
    int i, j, w, result;
    slong tmul = 10;

    /* Check f*g/g = f */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h, k, r;
       slong len, len1, len2;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;
       fmpz_mpoly_struct * qarr[1], * darr[1];

       fmpz_mpoly_ctx_init_rand(ctx, state, 10);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);
       fmpz_mpoly_init(k, ctx);
       fmpz_mpoly_init(r, ctx);

       len = n_randint(state, 100);
       len1 = n_randint(state, 100);
       len2 = n_randint(state, 100) + 1;

       exp_bits = n_randint(state, 200) + 1;
       exp_bits1 = n_randint(state, 200) + 1;
       exp_bits2 = n_randint(state, 200) + 1;

       coeff_bits = n_randint(state, 200);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
          do {
             fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits + 1, exp_bits2, ctx);
          } while (g->length == 0);
          fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
          fmpz_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);
          fmpz_mpoly_randtest_bits(r, state, len, coeff_bits, exp_bits, ctx);

          fmpz_mpoly_mul_johnson(h, f, g, ctx);

          qarr[0] = k;
          darr[0] = g;

          fmpz_mpoly_divrem_ideal_monagan_pearce(qarr, r, h, darr, 1, ctx);

          result = fmpz_mpoly_equal(f, k, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(k, ctx);
       fmpz_mpoly_clear(r, ctx);
    }

    /* Check f = g1*q1 + ... + gn*qn + r for random polys */
    for (i = 0; i < tmul*tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, r, k1, k2;
        fmpz_mpoly_struct * g, * q;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
        slong coeff_bits;
        slong n;
        fmpz_mpoly_struct * qarr[5], * darr[5];
        fmpz * shifts, * strides;

        num = n_randint(state, 5) + 1;

        g = (fmpz_mpoly_struct *) flint_malloc(num*sizeof(fmpz_mpoly_struct));
        q = (fmpz_mpoly_struct *) flint_malloc(num*sizeof(fmpz_mpoly_struct));

        fmpz_mpoly_ctx_init_rand(ctx, state, 9);

        for (w = 0; w < num; w++)
        {
            fmpz_mpoly_init(g + w, ctx);
            darr[w] = g + w;

            fmpz_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(k1, ctx);
        fmpz_mpoly_init(k2, ctx);
        fmpz_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound =  n_randint(state, 2 + 175/n/n) + 1;
        exp_bound1 = n_randint(state, 2 + 175/n/n) + 1;
        exp_bound2 = n_randint(state, 2 + 175/n/n) + 1;

        coeff_bits = n_randint(state, 70);

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
            fmpz_mpoly_inflate(f, f, shifts, strides, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);
            for (w = 0; w < num; w++)
            {
                do {
                    fmpz_mpoly_randtest_bound(darr[w], state, len2, coeff_bits + 1, exp_bound2, ctx);
                } while (darr[w]->length == 0);
                fmpz_mpoly_inflate(darr[w], darr[w], shifts, strides, ctx);
                fmpz_mpoly_assert_canonical(darr[w], ctx);
                fmpz_mpoly_randtest_bound(qarr[w], state, len, coeff_bits, exp_bound, ctx);
            }
            fmpz_mpoly_randtest_bound(k1, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bound(k2, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_divrem_ideal_monagan_pearce(qarr, r, f, darr, num, ctx);

            fmpz_mpoly_zero(k2, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mpoly_assert_canonical(qarr[w], ctx);
                fmpz_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
                fmpz_mpoly_add(k2, k2, k1, ctx);
            }
            fmpz_mpoly_add(k2, k2, r, ctx);

            result = fmpz_mpoly_equal(f, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = g1*q1 + ... + gn*qn + r for random polys"
                                               "\ni = %wd, j = %wd\n", i, j);
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

        for (w = 0; w < num; w++)
            fmpz_mpoly_clear(qarr[w], ctx);
        for (w = 0; w < num; w++)
            fmpz_mpoly_clear(darr[w], ctx);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(k1, ctx);
        fmpz_mpoly_clear(k2, ctx);
        fmpz_mpoly_clear(r, ctx);

        flint_free(g);
        flint_free(q);
    }

    /* Check aliasing */
    for (i = 0; i < 2*tmul * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, r, k1, k2;
       fmpz_mpoly_struct * g, * q;
       slong len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
       slong coeff_bits;
       slong n;
       fmpz_mpoly_struct * qarr[5], * darr[5];

       num = n_randint(state, 5) + 1;

       g = (fmpz_mpoly_struct *) flint_malloc(num*sizeof(fmpz_mpoly_struct));
       q = (fmpz_mpoly_struct *) flint_malloc(num*sizeof(fmpz_mpoly_struct));

       fmpz_mpoly_ctx_init_rand(ctx, state, 10);

       for (w = 0; w < num; w++)
       {
          fmpz_mpoly_init(g + w, ctx);
          darr[w] = g + w;

          fmpz_mpoly_init(q + w, ctx);
          qarr[w] = q + w;
       }

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(k1, ctx);
       fmpz_mpoly_init(k2, ctx);
       fmpz_mpoly_init(r, ctx);

       len = n_randint(state, 10);
       len1 = n_randint(state, 10);
       len2 = n_randint(state, 10) + 1;

       n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
       exp_bound =  n_randint(state, 3 + 200/n/n) + 1;
       exp_bound1 = n_randint(state, 3 + 200/n/n) + 1;
       exp_bound2 = n_randint(state, 3 + 200/n/n) + 1;

       coeff_bits = n_randint(state, 70);

       for (j = 0; j < 4; j++)
       {
          fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
          for (w = 0; w < num; w++)
          {
             do {
                fmpz_mpoly_randtest_bound(darr[w], state, len2, coeff_bits + 1, exp_bound2, ctx);
             } while (darr[w]->length == 0);
             fmpz_mpoly_randtest_bound(qarr[w], state, len, coeff_bits, exp_bound, ctx);
          }
          fmpz_mpoly_randtest_bound(k1, state, len, coeff_bits, exp_bound, ctx);
          fmpz_mpoly_randtest_bound(k2, state, len, coeff_bits, exp_bound, ctx);

          fmpz_mpoly_set(r, f, ctx);

          fmpz_mpoly_divrem_ideal_monagan_pearce(qarr, f, f, darr, num, ctx);

          fmpz_mpoly_zero(k2, ctx);
          for (w = 0; w < num; w++)
          {
             fmpz_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
             fmpz_mpoly_add(k2, k2, k1, ctx);
	      }
          fmpz_mpoly_add(k2, k2, f, ctx);

          result = fmpz_mpoly_equal(r, k2, ctx);

          if (!result)
          {
             printf("FAIL\n");
             flint_printf("Check aliasing\ni = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       for (w = 0; w < num; w++)
          fmpz_mpoly_clear(qarr[w], ctx);
       for (w = 0; w < num; w++)
          fmpz_mpoly_clear(darr[w], ctx);
       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(k1, ctx);
       fmpz_mpoly_clear(k2, ctx);
       fmpz_mpoly_clear(r, ctx);

       flint_free(g);
       flint_free(q);
    }

    TEST_FUNCTION_END(state);
}
