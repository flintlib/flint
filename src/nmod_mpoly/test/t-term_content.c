/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_term_content, state)
{
    int i, j, result;

    /* Check division by content leaves trivial content */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, k;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        nmod_mpoly_ctx_init_rand(ctx, state, 20, 2 + n_randint(state, -UWORD(2)));

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50) + 1;

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        nmod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);

        for (j = 0; j < 10; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);

            nmod_mpoly_term_content(f, g, ctx);

            if (nmod_mpoly_is_zero(g, ctx))
            {
                result = nmod_mpoly_is_zero(f, ctx);
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
                result = nmod_mpoly_length(f, ctx) == 1;
                if (!result)
                {
                    printf("FAIL\n");
                    flint_printf("Check content is monomial\ni = %wd, j = %wd\n", i ,j);
                    fflush(stdout);
                    flint_abort();
                }

                result = f->coeffs[0] == 1;
                if (!result)
                {
                    printf("FAIL\n");
                    flint_printf("Check content is monic\ni = %wd, j = %wd\n", i ,j);
                    fflush(stdout);
                    flint_abort();
                }

                result = nmod_mpoly_divides(k, g, f, ctx);
                if (!result)
                {
                    printf("FAIL\n");
                    flint_printf("Check content divides\ni = %wd, j = %wd\n", i ,j);
                    fflush(stdout);
                    flint_abort();
                }

                nmod_mpoly_term_content(k, k, ctx);
                result = nmod_mpoly_is_one(k, ctx);
                if (!result)
                {
                    printf("FAIL\n");
                    flint_printf("Check quotient is primitive\ni = %wd, j = %wd\n", i ,j);
                    fflush(stdout);
                    flint_abort();
                }
            }

            nmod_mpoly_term_content(g, g, ctx);
            result = nmod_mpoly_equal(f, g, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing \ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
