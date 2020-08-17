/*
    Copyright (C) 2018 Daniel Schultz

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
    int i, j, success;

    FLINT_TEST_INIT(state);

    flint_printf("inflate_deflate....");
    fflush(stdout);

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

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

