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

TEST_FUNCTION_START(fmpq_mpoly_get_set_is_fmpq, state)
{
    slong i;

    /* Set to random integer and compare */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f;
        fmpq_t c, d;
        slong len;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpq_init(c);
        fmpq_init(d);

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);

        len = n_randint(state, 10);
        exp_bits = n_randint(state, 200);
        coeff_bits = n_randint(state, 200);

        fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        if (fmpq_mpoly_is_fmpq(f, ctx))
        {
            fmpq_mpoly_get_fmpq(c, f, ctx);
            if (!fmpq_mpoly_equal_fmpq(f, c, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check is_fmpq and get_fmpq catch constants\ni = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_randtest(c, state, n_randint(state, 200) + 1);
        fmpq_mpoly_set_fmpq(f, c, ctx);
        if (!fmpq_mpoly_is_fmpq(f, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check set_fmpq makes is_fmpq true\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }
        fmpq_mpoly_get_fmpq(d, f, ctx);
        if (!fmpq_equal(c, d))
        {
            printf("FAIL\n");
            flint_printf("Check get_fmpq matches set_fmpq true\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpq_mpoly_clear(f, ctx);

        fmpq_clear(c);
        fmpq_clear(d);
    }

    TEST_FUNCTION_END(state);
}
