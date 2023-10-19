/*
    Copyright (C) 2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_neg, state)
{
    int i, result;

    /* Check -(-a) == a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       slong len, coeff_bits, exp_bits;

       fmpz_mpoly_ctx_init_rand(ctx, state, 20);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, 200) + 1;
       coeff_bits = n_randint(state, 200);

       fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

       fmpz_mpoly_neg(h, f, ctx);
       fmpz_mpoly_assert_canonical(h, ctx);
       fmpz_mpoly_neg(g, h, ctx);
       fmpz_mpoly_assert_canonical(g, ctx);

       result = fmpz_mpoly_equal(f, g, ctx);

       if (!result)
       {
             printf("FAIL\n");
             flint_printf("Check -(-a) == a\ni = %wd\n", i);
             fflush(stdout);
             flint_abort();
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
    }

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g;
       slong len, coeff_bits, exp_bits;

       fmpz_mpoly_ctx_init_rand(ctx, state, 20);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, 200) + 1;
       coeff_bits = n_randint(state, 200);

       fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

       fmpz_mpoly_neg(g, f, ctx);
       fmpz_mpoly_neg(g, g, ctx);

       result = fmpz_mpoly_equal(f, g, ctx);

       if (!result)
       {
             printf("FAIL\n");
             flint_printf("Check aliasing\ni = %wd\n", i);
             fflush(stdout);
             flint_abort();
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
    }

    TEST_FUNCTION_END(state);
}
