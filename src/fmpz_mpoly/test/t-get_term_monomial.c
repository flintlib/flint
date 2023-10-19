/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_get_term_monomial, state)
{
    int i, j;

    /* Check getting a coeff by its monomial */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_t c, d;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        flint_bitcnt_t exp_bits1, exp_bits2, exp_bits3;
        flint_bitcnt_t coeff_bits;
        slong len1, len2, len3;

        fmpz_init(c);
        fmpz_init(d);

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        len3 = n_randint(state, 100);

        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;
        exp_bits3 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
        fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
        fmpz_mpoly_randtest_bits(h, state, len3, coeff_bits, exp_bits3, ctx);

        fmpz_mpoly_repack_bits(h, f,
                                f->bits + n_randint(state, 2*FLINT_BITS), ctx);

        for (j = fmpz_mpoly_length(f, ctx) - 1; j >= 0; j--)
        {
            fmpz_mpoly_get_term_monomial(g, f, j, ctx);
            fmpz_mpoly_repack_bits(g, g,
                                  g->bits + n_randint(state, FLINT_BITS), ctx);
            fmpz_mpoly_get_term_coeff_fmpz(d, f, j, ctx);
            fmpz_mpoly_get_coeff_fmpz_monomial(c, h, g, ctx);

            if (!fmpz_equal(c, d))
            {
                flint_printf("FAIL\nCheck getting a coeff by its monomial\ni = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(c);
        fmpz_clear(d);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
