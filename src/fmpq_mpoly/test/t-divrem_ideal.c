/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_divrem_ideal, state)
{
    int i, j, w, result;

    /* Check s*f = g1*q1 + ... + gn*qn + r for random polys */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, r, k1, k2;
        fmpq_mpoly_struct * g, * q;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
        slong coeff_bits;
        fmpq_mpoly_struct * qarr[5], * darr[5];
        fmpz * shifts, * strides;

        num = n_randint(state, 5) + 1;

        g = (fmpq_mpoly_struct *) flint_malloc(num*sizeof(fmpq_mpoly_struct));
        q = (fmpq_mpoly_struct *) flint_malloc(num*sizeof(fmpq_mpoly_struct));

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = fmpq_mpoly_ctx_nvars(ctx);

        for (w = 0; w < num; w++)
        {
            fmpq_mpoly_init(g + w, ctx);
            darr[w] = g + w;
            fmpq_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(k1, ctx);
        fmpq_mpoly_init(k2, ctx);
        fmpq_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 12);
        len2 = n_randint(state, 8) + 1;

        exp_bound = n_randint(state, 10/FLINT_MAX(WORD(1), nvars) + 1) + 2;
        exp_bound1 = n_randint(state, 25/FLINT_MAX(WORD(1), nvars) + 1) + 2;
        exp_bound2 = n_randint(state, 20/FLINT_MAX(WORD(1), nvars) + 1) + 1;

        coeff_bits = n_randint(state, 40);

        shifts = (fmpz *) flint_malloc(nvars*sizeof(fmpz));
        strides = (fmpz *) flint_malloc(nvars*sizeof(fmpz));
        for (j = 0; j < nvars; j++)
        {
            fmpz_init(shifts + j);
            fmpz_init(strides + j);
            fmpz_randtest_unsigned(shifts + j, state, 100);
            fmpz_randtest_unsigned(strides + j, state, 100);
            fmpz_add_ui(strides + j, strides + j, 1);
        }

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpq_mpoly_inflate(f, f, shifts, strides, ctx);
            fmpq_mpoly_assert_canonical(f, ctx);
            for (w = 0; w < num; w++)
            {
                do {
                    fmpq_mpoly_randtest_bound(darr[w], state, len2, coeff_bits + 1, exp_bound2, ctx);
                } while (fmpq_mpoly_is_zero(darr[w], ctx));
                fmpq_mpoly_inflate(darr[w], darr[w], shifts, strides, ctx);
                fmpq_mpoly_assert_canonical(darr[w], ctx);
                fmpq_mpoly_randtest_bound(qarr[w], state, len, coeff_bits, exp_bound, ctx);
            }
            fmpq_mpoly_randtest_bound(k1, state, len, coeff_bits, exp_bound, ctx);
            fmpq_mpoly_randtest_bound(k2, state, len, coeff_bits, exp_bound, ctx);

            fmpq_mpoly_divrem_ideal(qarr, r, f, darr, num, ctx);
            fmpq_mpoly_assert_canonical(r, ctx);

            fmpq_mpoly_set(k2, r, ctx);
            for (w = 0; w < num; w++)
            {
                fmpq_mpoly_assert_canonical(qarr[w], ctx);
                fmpq_mpoly_remainder_test(r, darr[w], ctx);
                fmpq_mpoly_mul(k1, qarr[w], darr[w], ctx);
                fmpq_mpoly_add(k2, k2, k1, ctx);
	        }

            result = fmpq_mpoly_equal(f, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check s*f = g1*q1 + ... + gn*qn + r\ni=%wd j=%wd\n",i,j);
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

        for (w = 0; w < num; w++)
            fmpq_mpoly_clear(qarr[w], ctx);
        for (w = 0; w < num; w++)
            fmpq_mpoly_clear(darr[w], ctx);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(k1, ctx);
        fmpq_mpoly_clear(k2, ctx);
        fmpq_mpoly_clear(r, ctx);

        flint_free(g);
        flint_free(q);
    }

    /* Check aliasing of remainder */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, r, k1, k2;
        fmpq_mpoly_struct * g, * q;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
        slong coeff_bits;
        fmpq_mpoly_struct * qarr[5], * darr[5];

        num = n_randint(state, 5) + 1;

        g = (fmpq_mpoly_struct *) flint_malloc(num*sizeof(fmpq_mpoly_struct));
        q = (fmpq_mpoly_struct *) flint_malloc(num*sizeof(fmpq_mpoly_struct));

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = fmpq_mpoly_ctx_nvars(ctx);

        for (w = 0; w < num; w++)
        {
            fmpq_mpoly_init(g + w, ctx);
            darr[w] = g + w;
            fmpq_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(k1, ctx);
        fmpq_mpoly_init(k2, ctx);
        fmpq_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 12);
        len2 = n_randint(state, 8) + 1;

        exp_bound = n_randint(state, 10/FLINT_MAX(WORD(1), nvars) + 1) + 2;
        exp_bound1 = n_randint(state, 25/FLINT_MAX(WORD(1), nvars) + 1) + 2;
        exp_bound2 = n_randint(state, 20/FLINT_MAX(WORD(1), nvars) + 1) + 1;

        coeff_bits = n_randint(state, 40);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);

            for (w = 0; w < num; w++)
            {
                do {
                    fmpq_mpoly_randtest_bound(darr[w], state, len2, coeff_bits + 1, exp_bound2, ctx);
                } while (fmpq_mpoly_is_zero(darr[w], ctx));
                fmpq_mpoly_randtest_bound(qarr[w], state, len, coeff_bits, exp_bound, ctx);
            }
            fmpq_mpoly_randtest_bound(k1, state, len, coeff_bits, exp_bound, ctx);
            fmpq_mpoly_randtest_bound(k2, state, len, coeff_bits, exp_bound, ctx);

            fmpq_mpoly_set(r, f, ctx);
            fmpq_mpoly_divrem_ideal(qarr, r, r, darr, num, ctx);
            fmpq_mpoly_assert_canonical(r, ctx);

            fmpq_mpoly_set(k2, r, ctx);
            for (w = 0; w < num; w++)
            {
                fmpq_mpoly_assert_canonical(qarr[w], ctx);
                fmpq_mpoly_remainder_test(r, darr[w], ctx);
                fmpq_mpoly_mul(k1, qarr[w], darr[w], ctx);
                fmpq_mpoly_add(k2, k2, k1, ctx);
	        }

            result = fmpq_mpoly_equal(f, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of remainder\ni=%wd j=%wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (w = 0; w < num; w++)
            fmpq_mpoly_clear(qarr[w], ctx);
        for (w = 0; w < num; w++)
            fmpq_mpoly_clear(darr[w], ctx);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(k1, ctx);
        fmpq_mpoly_clear(k2, ctx);
        fmpq_mpoly_clear(r, ctx);

        flint_free(g);
        flint_free(q);
    }

    TEST_FUNCTION_END(state);
}
