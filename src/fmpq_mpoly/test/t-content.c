/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_content, state)
{
    slong i, j, k;

    /* Check content is gcd of coefficients */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        fmpq_t g1, g2, t;
        slong len;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);
        fmpq_init(g1);
        fmpq_init(g2);
        fmpq_init(t);

        len = n_randint(state, 50);
        exp_bits = n_randint(state, 200) + 2;
        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 10; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

            fmpq_zero(g1);
            for (k = 0; k < fmpq_mpoly_length(f, ctx); k++)
            {
                fmpq_mpoly_get_term_coeff_fmpq(t, f, k, ctx);
                fmpq_gcd(g1, g1, t);
            }

            fmpq_mpoly_content(g2, f, ctx);

            if (!fmpq_equal(g1, g2))
            {
                printf("FAIL\n");
                flint_printf("Check content\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_clear(g1);
        fmpq_clear(g2);
        fmpq_clear(t);
        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
