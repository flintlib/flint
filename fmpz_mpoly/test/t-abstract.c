/*
    Copyright (C) 2017 William Hart

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
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("abstract....");
    fflush(stdout);

    /* Check f = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, r;
        const fmpz_mpoly_struct * fptr[] = {f};
        slong len;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(r, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 2;
        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 1; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_abstract_add(r, (void **) fptr, 1, f->bits, ctx, NULL, NULL);
            result = fmpz_mpoly_equal(f, r, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = f\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(r, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check f + g = g + f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        const fmpz_mpoly_struct * fptr[] = {f, g};
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_add(h, f, g, ctx);
            fmpz_mpoly_abstract_add(k, (void **) fptr, 2, FLINT_MAX(f->bits, g->bits), ctx, NULL, NULL);
            result = fmpz_mpoly_equal(h, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f + g = g + f\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check f + (g + h) = (f + g) + h */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k1, k2;
        const fmpz_mpoly_struct * fptr[] = {f, g, h};
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k1, ctx);
        fmpz_mpoly_init(k2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(k1, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(k2, state, len, coeff_bits, exp_bits, ctx);

            fmpz_mpoly_add(k1, f, g, ctx);
            fmpz_mpoly_add(k1, k1, h, ctx);
            fmpz_mpoly_abstract_add(k2, (void **) fptr, 3, FLINT_MAX(FLINT_MAX(f->bits, g->bits), h->bits), ctx, NULL, NULL);
            result = fmpz_mpoly_equal(k1, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f + (g + h) = (f + g) + h\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k1, ctx);
        fmpz_mpoly_clear(k2, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
