/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_get_set_str_pretty, state)
{
    slong i;

    {
        slong len1;
        flint_bitcnt_t exp_bits1;
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, f1;
        char * str;
        const char * vars[] = {"x","y","z","w","u","v"};

        for (i = 0; i < 10*flint_test_multiplier(); i++)
        {
            fq_nmod_mpoly_ctx_init_rand(ctx, state, 6, FLINT_BITS, 10);

            fq_nmod_mpoly_init(f, ctx);
            fq_nmod_mpoly_init(f1, ctx);

            for (len1 = 3; len1 < 100; len1 += len1/2)
            {
                exp_bits1 = 3 + n_randint(state, 200);
                fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
                fq_nmod_mpoly_assert_canonical(f, ctx);

                str = fq_nmod_mpoly_get_str_pretty(f, vars, ctx);
                fq_nmod_mpoly_set_str_pretty(f1, str, vars, ctx);
                flint_free(str);

                if (!fq_nmod_mpoly_equal(f, f1, ctx))
                {
                    flint_printf("FAIL\n");
                    fflush(stdout);
                    flint_abort();
                }

            }

            fq_nmod_mpoly_clear(f, ctx);
            fq_nmod_mpoly_clear(f1, ctx);

            fq_nmod_mpoly_ctx_clear(ctx);
        }
    }

    TEST_FUNCTION_END(state);
}
