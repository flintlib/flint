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
#include "fmpq_mpoly.h"

int
main(void)
{
    slong i, j, tmul = 10;
    FLINT_TEST_INIT(state);

    flint_printf("div....");
    fflush(stdout);

    {
        fmpq_mpoly_t f, g, p, q;
        fmpq_mpoly_ctx_t ctx;
        const char * vars[] = {"x", "y", "z", "t", "u"};

        fmpq_mpoly_ctx_init(ctx, 5, ORD_LEX);
        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(p, ctx);
        fmpq_mpoly_init(q, ctx);

        fmpq_mpoly_set_str_pretty(f, "(1+x+y+2*z^2+3*t^3+5*u^5)^6", vars, ctx);
        fmpq_mpoly_set_str_pretty(g, "(1+u+t+2*z^2+3*y^3+5*x^5)^6", vars, ctx);

        fmpq_mpoly_mul(p, f, g, ctx);
        fmpq_mpoly_assert_canonical(p, ctx);

        fmpq_mpoly_div(q, p, g, ctx);
        fmpq_mpoly_assert_canonical(q, ctx);

        if (!fmpq_mpoly_equal(q, f, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check example\n");
            flint_abort();
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(p, ctx);
        fmpq_mpoly_clear(q, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f */
    for (i = 0; i < 10 * tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        flint_bitcnt_t coeff_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50) + 1;

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

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
            fmpq_mpoly_div(k, h, g, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);

            if (!fmpq_mpoly_equal(k, f, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni=%wd j=%wd\n",i,j);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(k, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check output agrees with divrem for random polys */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k, r;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        flint_bitcnt_t coeff_bits;
        fmpz * shifts, * strides;
        slong nvars;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = ctx->zctx->minfo->nvars;

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);
        fmpq_mpoly_init(r, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        coeff_bits = n_randint(state, 70);

        max_bound = 1 + 400/nvars/nvars;
        exp_bound = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        shifts = (fmpz *) flint_malloc(nvars*sizeof(fmpz));
        strides = (fmpz *) flint_malloc(nvars*sizeof(fmpz));
        for (j = 0; j < nvars; j++)
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
            fmpq_mpoly_randtest_bounds(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpq_mpoly_randtest_bounds(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));

            fmpq_mpoly_randtest_bounds(h, state, len, coeff_bits, exp_bound, ctx);
            fmpq_mpoly_randtest_bounds(k, state, len, coeff_bits, exp_bound, ctx);

            fmpq_mpoly_inflate(f, f, shifts, strides, ctx);
            fmpq_mpoly_inflate(g, g, shifts, strides, ctx);

            fmpq_mpoly_divrem(h, r, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_div(k, f, g, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);

            if (!fmpq_mpoly_equal(k, h, ctx))
            {
                printf("FAIL\n");                   
                flint_printf("Check output agrees with divrem\n"
                                                   "i = %wd, j = %wd\n", i, j);
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

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(k, ctx);
        fmpq_mpoly_clear(r, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing of quotient with first argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, r;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        flint_bitcnt_t coeff_bits;
        fmpz * shifts, * strides;
        slong nvars;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = ctx->zctx->minfo->nvars;

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(r, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        coeff_bits = n_randint(state, 70);

        max_bound = 1 + 400/nvars/nvars;
        exp_bound = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        shifts = (fmpz *) flint_malloc(nvars*sizeof(fmpz));
        strides = (fmpz *) flint_malloc(nvars*sizeof(fmpz));
        for (j = 0; j < nvars; j++)
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
            fmpq_mpoly_randtest_bounds(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpq_mpoly_randtest_bounds(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));

            fmpq_mpoly_randtest_bounds(h, state, len, coeff_bits, exp_bound, ctx);

            fmpq_mpoly_inflate(f, f, shifts, strides, ctx);
            fmpq_mpoly_inflate(g, g, shifts, strides, ctx);

            fmpq_mpoly_divrem(h, r, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_div(f, f, g, ctx);
            fmpq_mpoly_assert_canonical(f, ctx);

            if (!fmpq_mpoly_equal(f, h, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of quotient with first argument\n"
                                                   "i = %wd, j = %wd\n", i, j);
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

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(r, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing of quotient with second argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, r;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        flint_bitcnt_t coeff_bits;
        fmpz * shifts, * strides;
        slong nvars;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);
        nvars = ctx->zctx->minfo->nvars;

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(r, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        coeff_bits = n_randint(state, 70);

        max_bound = 1 + 400/nvars/nvars;
        exp_bound = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(nvars*sizeof(mp_limb_t));
        shifts = (fmpz *) flint_malloc(nvars*sizeof(fmpz));
        strides = (fmpz *) flint_malloc(nvars*sizeof(fmpz));
        for (j = 0; j < nvars; j++)
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
            fmpq_mpoly_randtest_bounds(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpq_mpoly_randtest_bounds(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));

            fmpq_mpoly_randtest_bounds(h, state, len, coeff_bits, exp_bound, ctx);

            fmpq_mpoly_inflate(f, f, shifts, strides, ctx);
            fmpq_mpoly_inflate(g, g, shifts, strides, ctx);

            fmpq_mpoly_divrem(h, r, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_div(f, f, g, ctx);
            fmpq_mpoly_assert_canonical(f, ctx);

            if (!fmpq_mpoly_equal(f, h, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of quotient with second argument\n"
                                                   "i = %wd, j = %wd\n", i, j);
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

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(r, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

