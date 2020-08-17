/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mpoly.h"

int
main(void)
{
    int i, j;
    FLINT_TEST_INIT(state);

    flint_printf("get_term....");
    fflush(stdout);

    /* Check a polynomial is the sum of its terms */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2, exp_bits3;
        slong len1, len2, len3;

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

        fmpz_mpoly_zero(h, ctx);
        for (j = fmpz_mpoly_length(f, ctx) - 1; j >= 0; j--)
        {
            fmpz_mpoly_get_term(g, f, j, ctx);
            fmpz_mpoly_add(h, h, g, ctx);
        }

        if (!fmpz_mpoly_equal(f, h, ctx))
        {
            flint_printf("FAIL\nCheck a polynomial is the sum of its terms\ni = %wd\n", i);
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

