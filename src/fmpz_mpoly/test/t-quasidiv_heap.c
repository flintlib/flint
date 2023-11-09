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

TEST_FUNCTION_START(fmpz_mpoly_quasidiv_heap, state)
{
    int result;
    slong i, j, tmul = 10;

    /* Check f*g/g = f */
    for (i = 0; i < 10 * tmul * flint_test_multiplier(); i++)
    {
        fmpz_t s1;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        flint_bitcnt_t coeff_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_init(s1);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100) + 1;

        exp_bits =  n_randint(state, 200) + 1;
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

            fmpz_mpoly_mul(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);

            fmpz_mpoly_quasidiv_heap(s1, k, h, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);

            result = fmpz_equal_ui(s1, 1) && fmpz_mpoly_equal(k, f, ctx);

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
        fmpz_clear(s1);
    }

    /* Check output agrees with divrem for random polys */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_t s1, s2;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k, r;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        flint_bitcnt_t coeff_bits;
        fmpz * shifts, * strides;
        slong n;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_init(s1);
        fmpz_init(s2);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);
        fmpz_mpoly_init(r, ctx);

        len = n_randint(state, 16);
        len1 = n_randint(state, 16);
        len2 = n_randint(state, 16) + 1;

        coeff_bits = n_randint(state, 70);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 200/n/n;
        exp_bound = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        shifts = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        strides = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            exp_bound[j] = UWORD(1) << (FLINT_BITS - 1);
            exp_bound1[j] = n_randint(state, max_bound) + 1;
            exp_bound2[j] = n_randint(state, max_bound) + 1;
            fmpz_init(shifts + j);
            fmpz_init(strides + j);
            fmpz_randtest_unsigned(shifts + j, state, 100);
            fmpz_randtest_unsigned(strides + j, state, 100);
            fmpz_add_ui(strides + j, strides + j, 1);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bounds(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpz_mpoly_randtest_bounds(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bounds(h, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bounds(k, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_inflate(f, f, shifts, strides, ctx);
            fmpz_mpoly_inflate(g, g, shifts, strides, ctx);

            fmpz_mpoly_quasidivrem_heap(s1, h, r, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_quasidiv_heap(s2, k, f, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);

            result = fmpz_equal(s1, s2) && fmpz_mpoly_equal(k, h, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check output agrees with divrem\n"
                                                   "i = %wd, j = %wd\n", i, j);
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

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_clear(r, ctx);
        fmpz_clear(s1);
        fmpz_clear(s2);
    }

    /* Check aliasing of quotient with first argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_t s1, s2;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, r;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        flint_bitcnt_t coeff_bits;
        fmpz * shifts, * strides;
        slong n;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_init(s1);
        fmpz_init(s2);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(r, ctx);

        len = n_randint(state, 16);
        len1 = n_randint(state, 16);
        len2 = n_randint(state, 16) + 1;

        coeff_bits = n_randint(state, 70);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 200/n/n;
        exp_bound = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        shifts = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        strides = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            exp_bound[j] = UWORD(1) << (FLINT_BITS - 1);
            exp_bound1[j] = n_randint(state, max_bound) + 1;
            exp_bound2[j] = n_randint(state, max_bound) + 1;
            fmpz_init(shifts + j);
            fmpz_init(strides + j);
            fmpz_randtest_unsigned(shifts + j, state, 100);
            fmpz_randtest_unsigned(strides + j, state, 100);
            fmpz_add_ui(strides + j, strides + j, 1);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bounds(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpz_mpoly_randtest_bounds(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (g->length == 0);

            fmpz_mpoly_randtest_bounds(h, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_inflate(f, f, shifts, strides, ctx);
            fmpz_mpoly_inflate(g, g, shifts, strides, ctx);

            fmpz_mpoly_quasidivrem_heap(s1, h, r, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_quasidiv_heap(s2, f, f, g, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);

            result = fmpz_equal(s1, s2) && fmpz_mpoly_equal(f, h, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of quotient with first argument\n"
                                                   "i = %wd, j = %wd\n", i, j);
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

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(r, ctx);
        fmpz_clear(s1);
        fmpz_clear(s2);
    }

    /* Check aliasing of quotient with second argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_t s1, s2;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, r;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        flint_bitcnt_t coeff_bits;
        fmpz * shifts, * strides;
        slong n;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_init(s1);
        fmpz_init(s2);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(r, ctx);

        len = n_randint(state, 16);
        len1 = n_randint(state, 16);
        len2 = n_randint(state, 16) + 1;

        coeff_bits = n_randint(state, 70);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 200/n/n;
        exp_bound = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        shifts = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        strides = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            exp_bound[j] = UWORD(1) << (FLINT_BITS - 1);
            exp_bound1[j] = n_randint(state, max_bound) + 1;
            exp_bound2[j] = n_randint(state, max_bound) + 1;
            fmpz_init(shifts + j);
            fmpz_init(strides + j);
            fmpz_randtest_unsigned(shifts + j, state, 100);
            fmpz_randtest_unsigned(strides + j, state, 100);
            fmpz_add_ui(strides + j, strides + j, 1);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bounds(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpz_mpoly_randtest_bounds(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (g->length == 0);

            fmpz_mpoly_randtest_bounds(h, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_inflate(f, f, shifts, strides, ctx);
            fmpz_mpoly_inflate(g, g, shifts, strides, ctx);

            fmpz_mpoly_quasidivrem_heap(s1, h, r, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_quasidiv_heap(s2, g, f, g, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);

            result = fmpz_equal(s1, s2) && fmpz_mpoly_equal(g, h, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of quotient with second argument\n"
                                                   "i = %wd, j = %wd\n", i, j);
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

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(r, ctx);
        fmpz_clear(s1);
        fmpz_clear(s2);
    }

    TEST_FUNCTION_END(state);
}
