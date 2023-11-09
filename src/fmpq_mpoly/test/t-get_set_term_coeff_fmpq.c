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

TEST_FUNCTION_START(fmpq_mpoly_get_set_term_coeff_fmpq, state)
{
    slong i, j;

    /* check get and set match */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        fmpq_t c, d;
        slong len, coeff_bits, exp_bits, index;

        fmpq_init(c);
        fmpq_init(d);
        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);

        do {
            len = n_randint(state, 100);
            exp_bits = n_randint(state, 200) + 1;
            coeff_bits = n_randint(state, 200);

            fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
        } while (fmpq_mpoly_length(f, ctx) == 0);

        for (j = 0; j < 10; j++)
        {
            fmpq_randtest(c, state, n_randint(state, 100) + 1);

            index = n_randint(state, fmpq_mpoly_length(f, ctx));
            fmpq_mpoly_set_term_coeff_fmpq(f, index, c, ctx);
            fmpq_mpoly_get_term_coeff_fmpq(d, f, index, ctx);
            if (!fmpq_equal(c, d))
            {
                printf("FAIL\n");
                flint_printf("check get and set match\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpq_mul_fmpz(c, fmpq_mpoly_content_ref(f, ctx),
                        fmpq_mpoly_zpoly_term_coeff_ref(f, index, ctx));
            if (!fmpq_equal(c, d))
            {
                printf("FAIL\n");
                flint_printf("check reference match\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_clear(c);
        fmpq_clear(d);
    }

    TEST_FUNCTION_END(state);
}
