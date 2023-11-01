/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_scalar_mul_fmpz, state)
{
    slong i, j;

    /* Check (f*a)*b = f*(a*b) */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
       fmpz_mod_mpoly_ctx_t ctx;
       fmpz_mod_mpoly_t f, g, h, k;
       fmpz_t a, b, c;
       slong len, exp_bits;

       fmpz_init(a);
       fmpz_init(b);
       fmpz_init(c);

       fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);

       fmpz_mod_mpoly_init(f, ctx);
       fmpz_mod_mpoly_init(g, ctx);
       fmpz_mod_mpoly_init(h, ctx);
       fmpz_mod_mpoly_init(k, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, 200) + 1;

       for (j = 0; j < 10; j++)
       {
          fmpz_mod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
          fmpz_mod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);
          fmpz_mod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
          fmpz_mod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);

          fmpz_randtest(a, state, n_randint(state, 200));
          fmpz_randtest(b, state, n_randint(state, 200));
          fmpz_mul(c, a, b);

          fmpz_mod_mpoly_scalar_mul_fmpz(g, f, a, ctx);
          fmpz_mod_mpoly_scalar_mul_fmpz(h, g, b, ctx);

          fmpz_mod_mpoly_scalar_mul_fmpz(k, f, c, ctx);

          if (!fmpz_mod_mpoly_equal(h, k, ctx))
          {
             flint_printf("FAIL: Check (f*a)*b = f*(a*b)\n");
             flint_printf("i = %wd, j = %wd\n", i, j);
             fflush(stdout);
             flint_abort();
          }
       }

       fmpz_mod_mpoly_clear(f, ctx);
       fmpz_mod_mpoly_clear(g, ctx);
       fmpz_mod_mpoly_clear(h, ctx);
       fmpz_mod_mpoly_clear(k, ctx);
       fmpz_mod_mpoly_ctx_clear(ctx);

       fmpz_clear(a);
       fmpz_clear(b);
       fmpz_clear(c);
    }

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h;
        fmpz_t c;
        slong len, exp_bits;

        fmpz_init(c);

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);

        len = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;

        for (j = 0; j < 10; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            fmpz_randtest(c, state, n_randint(state, 200));

            fmpz_mod_mpoly_set(g, f, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz(h, f, c, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz(g, g, c, ctx);

            if (!fmpz_mod_mpoly_equal(g, h, ctx))
            {
                flint_printf("FAIL: Check aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
       }

       fmpz_mod_mpoly_clear(f, ctx);
       fmpz_mod_mpoly_clear(g, ctx);
       fmpz_mod_mpoly_clear(h, ctx);
       fmpz_mod_mpoly_ctx_clear(ctx);

       fmpz_clear(c);
    }

    TEST_FUNCTION_END(state);
}
