/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_inflate_deflate, state)
{
    int i, j, success;

    /* Check deflate undoes inflate */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        slong k;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        fmpz * strides, * shifts;
        slong len1, len2, len3;
        flint_bitcnt_t exp_bits1, exp_bits2, exp_bits3;
        flint_bitcnt_t coeff_bits, stride_bits, shift_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);
        len3 = n_randint(state, 50);
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;
        exp_bits3 = n_randint(state, 100) + 2;
        coeff_bits = n_randint(state, 100);

        stride_bits = n_randint(state, 100) + 2;
        shift_bits = n_randint(state, 100) + 2;

        strides = flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        shifts = flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_init(strides + k);
            fmpz_init(shifts + k);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len3, coeff_bits, exp_bits3, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts + k, state, shift_bits);
                fmpz_randtest_not_zero(strides + k, state, stride_bits);
                fmpz_abs(strides + k, strides + k);
            }

            fmpz_mpoly_inflate(h, f, shifts, strides, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_deflate(g, h, shifts, strides, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

            if (!fmpz_mpoly_equal(f, g, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check deflate undoes inflate\n"
                                                     "i: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mpoly_set(h, f, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_inflate(h, h, shifts, strides, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_deflate(h, h, shifts, strides, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);

            if (!fmpz_mpoly_equal(f, h, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check deflate undoes inflate with aliasing\n"
                                                     "i: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_clear(strides + k);
            fmpz_clear(shifts + k);
        }
        flint_free(strides);
        flint_free(shifts);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check deflating by deflation leaves trivial deflation */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        slong k;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        fmpz * strides, * shifts;
        slong len1, len2, len3;
        flint_bitcnt_t exp_bits1, exp_bits2, exp_bits3;
        flint_bitcnt_t coeff_bits, stride_bits, shift_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);
        len3 = n_randint(state, 50);
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;
        exp_bits3 = n_randint(state, 100) + 2;
        coeff_bits = n_randint(state, 100);

        stride_bits = n_randint(state, 10) + 2;
        shift_bits = n_randint(state, 10) + 2;

        strides = flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        shifts = flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_init(strides + k);
            fmpz_init(shifts + k);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len3, coeff_bits, exp_bits3, ctx);

            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                fmpz_randtest_unsigned(shifts + k, state, shift_bits);
                fmpz_randtest_unsigned(strides + k, state, stride_bits);
            }

            fmpz_mpoly_inflate(h, f, shifts, strides, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_deflation(shifts, strides, h, ctx);
            fmpz_mpoly_deflate(g, h, shifts, strides, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);
            fmpz_mpoly_deflation(shifts, strides, g, ctx);
            fmpz_mpoly_deflate(f, g, shifts, strides, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);
            success = fmpz_mpoly_equal(f, g, ctx);
            for (k = 0; k < ctx->minfo->nvars; k++)
            {
                if (!fmpz_is_zero(shifts + k))
                    success = 0;
                if (fmpz_cmp_ui(strides + k, UWORD(1)) > 0)
                    success = 0;
            }
            if (!success)
            {
                printf("FAIL\n");
                flint_printf("Check deflating by deflation leaves trivial deflation\n"
                                                     "i: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (k = 0; k < ctx->minfo->nvars; k++)
        {
            fmpz_clear(strides + k);
            fmpz_clear(shifts + k);
        }
        flint_free(strides);
        flint_free(shifts);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check deflate with zero-stride does not divide by zero */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t x, res;
        fmpz *stride, *shift;

        fmpz_mpoly_ctx_init(ctx, 1, ORD_LEX);
        fmpz_mpoly_init(x, ctx);
        fmpz_mpoly_init(res, ctx);

        stride = flint_malloc(ctx->minfo->nvars * sizeof(fmpz));
        shift = flint_malloc(ctx->minfo->nvars * sizeof(fmpz));
        fmpz_init(stride + 0);
        fmpz_init(shift + 0);

        /* Set x to just the generator */
        fmpz_mpoly_gen(x, 0, ctx);

        /* --- With zero shift --- */
        fmpz_set_ui(shift + 0, 0);

        /* Attempt to deflate x to 1 with a shift and stride 0. */
        /* That is 1 -> (1 - 0) / 0 */
        /* Division by zero should not be raised here */
        fmpz_mpoly_deflate(res, x, shift, stride, ctx);

        if (!fmpz_mpoly_equal_ui(res, 1, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check deflate with zero-shift and zero-stride\n");
            fflush(stdout);
            flint_abort();
        }

        /* --- With non-zero shift --- */
        fmpz_set_ui(shift + 0, 1);

        /* Attempt to deflate x to 1 with a shift of 1 and stride 0. */
        /* That is 1 -> (1 - 1) / 0 */
        /* Division by zero should not be raised here */
        fmpz_mpoly_deflate(res, x, shift, stride, ctx);

        if (!fmpz_mpoly_equal_ui(res, 1, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check deflate with non-zero shift and zero-stride\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(stride + 0);
        fmpz_clear(shift + 0);
        flint_free(stride);
        flint_free(shift);

        fmpz_mpoly_clear(x, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
