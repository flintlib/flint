/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_quasidivrem_ideal_heap, state)
{
    int i, j, w, result;

    /* Check s*f = g1*q1 + ... + gn*qn + r for random polys */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t scale;
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

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        for (w = 0; w < num; w++)
        {
            fmpz_mpoly_init(g + w, ctx);
            darr[w] = g + w;
            fmpz_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }

        fmpz_init(scale);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(k1, ctx);
        fmpz_mpoly_init(k2, ctx);
        fmpz_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 12);
        len2 = n_randint(state, 8) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound = n_randint(state, 10/n + 1) + 2;
        exp_bound1 = n_randint(state, 25/n + 1) + 2;
        exp_bound2 = n_randint(state, 20/n + 1) + 1;

        coeff_bits = n_randint(state, 40);

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

            fmpz_mpoly_quasidivrem_ideal_heap(scale, qarr, r, f, darr, num, ctx);
            fmpz_mpoly_assert_canonical(r, ctx);

            fmpz_mpoly_set(k2, r, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mpoly_assert_canonical(qarr[w], ctx);
                fmpz_mpoly_remainder_strongtest(r, darr[w], ctx);
                fmpz_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
                fmpz_mpoly_add(k2, k2, k1, ctx);
	        }

            fmpz_mpoly_scalar_mul_fmpz(f, f, scale, ctx);
            result = fmpz_mpoly_equal(f, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check s*f = g1*q1 + ... + gn*qn + r\ni=%wd j=%wd\n",i,j);
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
        fmpz_clear(scale);
    }

    /* Check aliasing remainder */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t scale;
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

        fmpz_init(scale);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(k1, ctx);
        fmpz_mpoly_init(k2, ctx);
        fmpz_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 12);
        len2 = n_randint(state, 8) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound = n_randint(state, 10/n + 1) + 2;
        exp_bound1 = n_randint(state, 25/n + 1) + 2;
        exp_bound2 = n_randint(state, 20/n + 1) + 1;

        coeff_bits = n_randint(state, 40);

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
            fmpz_mpoly_quasidivrem_ideal_heap(scale, qarr, r, r, darr, num, ctx);
            fmpz_mpoly_assert_canonical(r, ctx);
            fmpz_mpoly_set(k2, r, ctx);
            for (w = 0; w < num; w++)
            {
                fmpz_mpoly_assert_canonical(qarr[w], ctx);
                fmpz_mpoly_remainder_strongtest(r, darr[w], ctx);
                fmpz_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
                fmpz_mpoly_add(k2, k2, k1, ctx);
	        }

            fmpz_mpoly_scalar_mul_fmpz(f, f, scale, ctx);
            result = fmpz_mpoly_equal(f, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing remainder\ni=%wd j=%wd\n",i,j);
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
        fmpz_clear(scale);
    }

    TEST_FUNCTION_END(state);
}
