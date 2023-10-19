/*
    Copyright (C) 2017 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_set_equal, state)
{
    int i, result;

    /* Set b = a and check a == b */
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

       fmpz_mpoly_set(g, f, ctx);

       result = fmpz_mpoly_equal(f, g, ctx);

       if (!result)
       {
          printf("FAIL\n");
          flint_printf("Set b = a and check a == b\ni = %wd\n", i);
          fflush(stdout);
          flint_abort();
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
    }

    /* Set b = a, alter b and check a == b */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g;
       ordering_t ord;
       slong nvars, len, coeff_bits, exp_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 20) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, 200) + 1;
       coeff_bits = n_randint(state, 200);

       fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

       fmpz_mpoly_add_ui(g, f, UWORD(1), ctx);

       result = !fmpz_mpoly_equal(f, g, ctx);

       if (!result)
       {
          printf("FAIL\n");
          flint_printf("Set b = a, alter b and check a == b\ni = %wd\n", i);
          fflush(stdout);
          flint_abort();
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
    }

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g;
       ordering_t ord;
       slong nvars, len, coeff_bits, exp_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 20) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, 200) + 1;
       coeff_bits = n_randint(state, 200);

       fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

       fmpz_mpoly_set(g, f, ctx);
       fmpz_mpoly_set(g, g, ctx);

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
