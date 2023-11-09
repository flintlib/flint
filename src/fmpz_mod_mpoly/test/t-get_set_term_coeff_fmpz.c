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

TEST_FUNCTION_START(fmpz_mod_mpoly_get_set_term_coeff_fmpz, state)
{
    slong i, j;

    /* Set coeff and get coeff and compare */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f;
        fmpz_t c, d;
        slong len, index;
        flint_bitcnt_t exp_bits;

        fmpz_init(c);
        fmpz_init(d);

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);
        fmpz_mod_mpoly_init(f, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 1;

        fmpz_mod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        if (f->length > 0)
        {
            for (j = 0; j < 10; j++)
            {
                fmpz_randtest(c, state, n_randint(state, 200));

                index = n_randint(state, f->length);

                fmpz_mod_mpoly_set_term_coeff_fmpz(f, index, c, ctx);
                fmpz_mod_mpoly_get_term_coeff_fmpz(d, f, index, ctx);
                if (!fmpz_mod_equal_fmpz(d, c, ctx->ffinfo))
                {
                    flint_printf("FAIL: check get and set match\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }

                if (!fmpz_equal(f->coeffs + index, d))
                {
                    flint_printf("FAIL: check reference match\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
