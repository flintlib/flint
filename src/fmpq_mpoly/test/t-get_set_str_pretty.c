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

TEST_FUNCTION_START(fmpq_mpoly_get_set_str_pretty, state)
{
    slong i;

    {
        slong len1, exp_bits, coeff_bits;
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, f1;
        char * str;
        const char * vars[] = {"x","xy","y","yx","z","zz"};

        for (i = 0; i < 1 * flint_test_multiplier(); i++)
        {
            fmpq_mpoly_ctx_init_rand(ctx, state, 6);
            fmpq_mpoly_init(f, ctx);
            fmpq_mpoly_init(f1, ctx);

            for (len1 = 3; len1 < 1000; len1 += len1/2)
            {
                coeff_bits = 100;
                exp_bits = 100;
                fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits, ctx);
                fmpq_mpoly_assert_canonical(f, ctx);
                str = fmpq_mpoly_get_str_pretty(f, vars, ctx);
                fmpq_mpoly_set_str_pretty(f1, str, vars, ctx);
                fmpq_mpoly_assert_canonical(f1, ctx);
                flint_free(str);

                if (!fmpq_mpoly_equal(f, f1, ctx))
                {
                    flint_printf("FAIL\n");
                    flint_printf("check that parsing inverts printing\ni = %wd, len1 = %wd\n", i ,len1);
                    fflush(stdout);
                    flint_abort();
                }
            }

            fmpq_mpoly_clear(f, ctx);
            fmpq_mpoly_clear(f1, ctx);
        }
    }

    TEST_FUNCTION_END(state);
}
