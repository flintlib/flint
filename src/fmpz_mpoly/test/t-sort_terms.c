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

TEST_FUNCTION_START(fmpz_mpoly_sort_terms, state)
{
    int i, j, result;

    /* Check scramble and sort */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g;
        slong len;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);

        len = n_randint(state, 200);
        exp_bits = n_randint(state, 200) + 1;
        coeff_bits = n_randint(state, 20);

        for (j = 0; j < 4; j++)
        {
            slong N, k;

            fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_set(g, f, ctx);

            N = mpoly_words_per_exp(f->bits, ctx->minfo);
            for (k = WORD(0); k < f->length; k++)
            {
                ulong a, b;
                a = n_randint(state, f->length);
                b = n_randint(state, f->length);
                fmpz_swap(f->coeffs + a, f->coeffs + b);
                mpoly_monomial_swap(f->exps + N*a, f->exps + N*b, N);
            }

            fmpz_mpoly_sort_terms(f, ctx);
            result = fmpz_mpoly_equal(f, g, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check scramble and sort\ni = %wd, j = %wd\n", i ,j);
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
