/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_repack_bits, state)
{
    int i, j, success;

    /* Check packing up */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g;
        slong len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2, newbits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;
        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);

            newbits = g->bits + n_randint(state, 2*FLINT_BITS);
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fmpz_mpoly_repack_bits(f, g, newbits, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);

            if (!success || !fmpz_mpoly_equal(f, g, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check packing up\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            newbits = g->bits + n_randint(state, 2*FLINT_BITS);
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fmpz_mpoly_repack_bits(g, g, newbits, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

            if (!success || !fmpz_mpoly_equal(f, g, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check packing up with aliasing\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check repacking down up */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2, newbits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;
        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len2, coeff_bits, exp_bits2, ctx);

            newbits = g->bits + n_randint(state, 2*FLINT_BITS);
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            fmpz_mpoly_repack_bits(f, g, newbits, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);

            newbits = g->bits + n_randint(state, 2*FLINT_BITS);
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fmpz_mpoly_repack_bits(h, f, newbits, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);

            if (!success || !fmpz_mpoly_equal(h, g, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check repacking down\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            newbits = g->bits + n_randint(state, 2*FLINT_BITS);
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fmpz_mpoly_repack_bits(f, f, newbits, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);

            if (!success || !fmpz_mpoly_equal(f, g, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check repacking down with aliasing\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check packing down */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g;
        slong len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2, newbits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);

        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;
        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);

            if (g->bits <= MPOLY_MIN_BITS)
                continue;

            newbits = n_randint(state, g->bits - MPOLY_MIN_BITS) + MPOLY_MIN_BITS;
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fmpz_mpoly_repack_bits(f, g, newbits, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);

            if (success && !fmpz_mpoly_equal(f, g, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check packing down\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mpoly_set(f, g, ctx);
            newbits = n_randint(state, g->bits - MPOLY_MIN_BITS) + MPOLY_MIN_BITS;
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fmpz_mpoly_repack_bits(g, g, newbits, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

            if (success && !fmpz_mpoly_equal(f, g, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check packing down with aliasing\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
