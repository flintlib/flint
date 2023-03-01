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
#include "fq_nmod_mpoly.h"


int
main(void)
{
    slong i, j;
    int success;

    FLINT_TEST_INIT(state);

    flint_printf("repack_bits....");
    fflush(stdout);

    /* Check packing up */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g;
        slong len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2, newbits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 9);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);

        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);

            newbits = g->bits + n_randint(state, 2*FLINT_BITS);
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fq_nmod_mpoly_repack_bits(f, g, newbits, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);

            if (!success || !fq_nmod_mpoly_equal(f, g, ctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check packing up\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            newbits = g->bits + n_randint(state, 2*FLINT_BITS);
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fq_nmod_mpoly_repack_bits(g, g, newbits, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);

            if (!success || !fq_nmod_mpoly_equal(f, g, ctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check packing up with aliasing\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check repacking down up */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        slong len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2, newbits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 9);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);

        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len2, exp_bits2, ctx);

            newbits = g->bits + n_randint(state, 2*FLINT_BITS);
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            fq_nmod_mpoly_repack_bits(f, g, newbits, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);

            newbits = g->bits + n_randint(state, 2*FLINT_BITS);
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fq_nmod_mpoly_repack_bits(h, f, newbits, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);

            if (!success || !fq_nmod_mpoly_equal(h, g, ctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check repacking down\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            newbits = g->bits + n_randint(state, 2*FLINT_BITS);
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fq_nmod_mpoly_repack_bits(f, f, newbits, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);

            if (!success || !fq_nmod_mpoly_equal(f, g, ctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check repacking down with aliasing\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check packing down */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g;
        slong len1, len2;
        flint_bitcnt_t exp_bits1, exp_bits2, newbits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 9);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);

        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);

            if (g->bits <= MPOLY_MIN_BITS)
                continue;

            newbits = n_randint(state, g->bits - MPOLY_MIN_BITS) + MPOLY_MIN_BITS;
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fq_nmod_mpoly_repack_bits(f, g, newbits, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);

            if (success && !fq_nmod_mpoly_equal(f, g, ctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check packing down\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_set(f, g, ctx);
            newbits = n_randint(state, g->bits - MPOLY_MIN_BITS) + MPOLY_MIN_BITS;
            newbits = mpoly_fix_bits(newbits, ctx->minfo);
            success = fq_nmod_mpoly_repack_bits(g, g, newbits, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);

            if (success && !fq_nmod_mpoly_equal(f, g, ctx))
            {
                flint_printf("FAIL\n");
                flint_printf("Check packing down with aliasing\ni: %wd  j: %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

