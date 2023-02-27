/*
    Copyright (C) 2017 William Hart
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
#include "ulong_extras.h"

int
main(void)
{
    slong i;
    FLINT_TEST_INIT(state);

    flint_printf("get_set_is_fmpz....");
    fflush(stdout);

    /* Set to random integer and compare */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f;
        fmpz_t c, d;
        slong len;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpz_init(c);
        fmpz_init(d);

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);

        len = n_randint(state, 10);
        exp_bits = n_randint(state, 200);
        coeff_bits = n_randint(state, 200);

        fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

        if (fmpz_mpoly_is_fmpz(f, ctx))
        {
            fmpz_mpoly_get_fmpz(c, f, ctx);
            if (!fmpz_mpoly_equal_fmpz(f, c, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check is_fmpz and get_fmpz catch constants\ni = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_randtest(c, state, n_randint(state, 200));
        fmpz_mpoly_set_fmpz(f, c, ctx);
        if (!fmpz_mpoly_is_fmpz(f, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check set_fmpz makes is_fmpz true\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }
        fmpz_mpoly_get_fmpz(d, f, ctx);
        if (!fmpz_equal(c, d))
        {
            printf("FAIL\n");
            flint_printf("Check get_fmpz matches set_fmpz true\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

