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

TEST_FUNCTION_START(fmpz_mpoly_term_content, state)
{
    int i, j, result;

    /* Check division by content leaves trivial content */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, k;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50) + 1;

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 20);

        fmpz_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);

            fmpz_mpoly_scalar_mul_ui(g, g, 1 + n_randint(state, 10), ctx);

            fmpz_mpoly_term_content(f, g, ctx);

            if (fmpz_mpoly_is_zero(g, ctx))
            {
                result = fmpz_mpoly_is_zero(f, ctx);
                if (!result)
                {
                    printf("FAIL\n");
                    flint_printf("Check zero\ni = %wd, j = %wd\n", i ,j);
                    fflush(stdout);
                    flint_abort();
                }
            }
            else
            {
                result = fmpz_mpoly_length(f, ctx) == 1;
                if (!result)
                {
                    printf("FAIL\n");
                    flint_printf("Check content is monomial\ni = %wd, j = %wd\n", i ,j);
                    fflush(stdout);
                    flint_abort();
                }

                result = fmpz_sgn(f->coeffs + 0) > 0;
                if (!result)
                {
                    printf("FAIL\n");
                    flint_printf("Check content has positive lc\ni = %wd, j = %wd\n", i ,j);
                    fflush(stdout);
                    flint_abort();
                }

                result = fmpz_mpoly_divides(k, g, f, ctx);
                if (!result)
                {
                    printf("FAIL\n");
                    flint_printf("Check content divides\ni = %wd, j = %wd\n", i ,j);
                    fflush(stdout);
                    flint_abort();
                }

                fmpz_mpoly_term_content(k, k, ctx);
                result = fmpz_mpoly_is_one(k, ctx);
                if (!result)
                {
                    printf("FAIL\n");
                    flint_printf("Check quotient is primitive\ni = %wd, j = %wd\n", i ,j);
                    fflush(stdout);
                    flint_abort();
                }
            }

            fmpz_mpoly_term_content(g, g, ctx);
            result = fmpz_mpoly_equal(f, g, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing \ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
