/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_divides, state)
{
    int i, j, result, ok1, ok2;

    /* Check f*g/g = f */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10) + 1;

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            do {
                fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits + 1, exp_bits2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

            fmpq_mpoly_mul(h, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            ok1 = fmpq_mpoly_divides(k, h, g, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);
            result = (ok1 && fmpq_mpoly_equal(f, k, ctx));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(k, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check random polys don't divide */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits, exp_bits1, exp_bits2, maxbits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        maxbits = 10/FLINT_MAX(WORD(1), fmpq_mpoly_ctx_nvars(ctx)) + 2;
        exp_bits = n_randint(state, maxbits) + 1;
        exp_bits1 = n_randint(state, maxbits) + 1;
        exp_bits2 = n_randint(state, maxbits) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpq_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));
            fmpq_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
            fmpq_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

            ok1 = fmpq_mpoly_divides(h, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);

            if (ok1)
            {
                fmpq_mpoly_mul(k, h, g, ctx);
                fmpq_mpoly_assert_canonical(k, ctx);
            }

            result = (ok1 == 0 || fmpq_mpoly_equal(f, k, ctx));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check random polys don't divide\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(k, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument, exact division */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50) + 1;

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            do {
                fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits + 1, exp_bits2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

            fmpq_mpoly_mul(h, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            ok1 = fmpq_mpoly_divides(k, h, g, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);
            ok2 = fmpq_mpoly_divides(h, h, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            result = (ok1 == 1 && ok2 == 1 && fmpq_mpoly_equal(h, k, ctx));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first argument, exact division\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(k, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing, first argument, random polys */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits, exp_bits1, exp_bits2, maxbits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        maxbits = 10/FLINT_MAX(WORD(1), fmpq_mpoly_ctx_nvars(ctx)) + 2;
        exp_bits = n_randint(state, maxbits) + 1;
        exp_bits1 = n_randint(state, maxbits) + 1;
        exp_bits2 = n_randint(state, maxbits) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpq_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));
            fmpq_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

            ok1 = fmpq_mpoly_divides(h, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            ok2 = fmpq_mpoly_divides(f, f, g, ctx);
            fmpq_mpoly_assert_canonical(f, ctx);

            result = ((ok1 == ok2) &&  (ok1 == 0 || fmpq_mpoly_equal(f, h, ctx)));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing, first argument, random polys\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument, exact division */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50) + 1;

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            do {
                fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits + 1, exp_bits2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

            fmpq_mpoly_mul(h, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            ok1 = fmpq_mpoly_divides(k, h, g, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);
            ok2 = fmpq_mpoly_divides(g, h, g, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);
            result = (ok1 == 1 && ok2 == 1 && fmpq_mpoly_equal(g, k, ctx));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second argument, exact division\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(k, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing, second argument, random polys */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits, exp_bits1, exp_bits2, maxbits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100) + 1;

        maxbits = 10/FLINT_MAX(WORD(1), fmpq_mpoly_ctx_nvars(ctx)) + 2;
        exp_bits = n_randint(state, maxbits) + 1;
        exp_bits1 = n_randint(state, maxbits) + 1;
        exp_bits2 = n_randint(state, maxbits) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpq_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));
            fmpq_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

            ok1 = fmpq_mpoly_divides(h, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            ok2 = fmpq_mpoly_divides(g, f, g, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);

            result = ((ok1 == ok2) &&  (ok1 == 0 || fmpq_mpoly_equal(g, h, ctx)));

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing, second argument, random polys\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
