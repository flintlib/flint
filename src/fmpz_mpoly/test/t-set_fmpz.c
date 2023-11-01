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

TEST_FUNCTION_START(fmpz_mpoly_set_fmpz, state)
{
    int i, result;

    /* Set to random integer and compare */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f;
       fmpz_t c, d;
       slong nvars, len, coeff_bits, exp_bits, j;
       ulong * exp;

       fmpz_init(c);
       fmpz_init(d);

       fmpz_mpoly_ctx_init_rand(ctx, state, 20);
       nvars = ctx->minfo->nvars;

       fmpz_mpoly_init(f, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, 200) + 1;
       coeff_bits = n_randint(state, 200);

       fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

       fmpz_randtest(c, state, n_randint(state, 200));

       fmpz_mpoly_set_fmpz(f, c, ctx);

       exp = (ulong *) flint_malloc(nvars*sizeof(ulong));

       for (j = 0; j < nvars; j++)
          exp[j] = 0;

       fmpz_mpoly_get_coeff_fmpz_ui(d, f, exp, ctx);

       result = fmpz_equal(c, d) &&
                ((fmpz_is_zero(c) && f->length == 0) || f->length == 1);

       if (!result)
       {
          printf("FAIL\n");
          flint_printf("Set to random integer and compare\ni = %wd\n", i);
          fflush(stdout);
          flint_abort();
       }

       flint_free(exp);

       fmpz_mpoly_clear(f, ctx);

       fmpz_clear(c);
       fmpz_clear(d);
    }

    TEST_FUNCTION_END(state);
}
