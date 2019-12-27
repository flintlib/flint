/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"

int
main(void)
{
    slong i, j;
    int success;

    FLINT_TEST_INIT(state);

    flint_printf("univar....");
    fflush(stdout);


    /* Check mpoly -> mpoly_univar -> mpoly */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g;
        fmpz_mpoly_univar_t fx;
        slong len1, len2;
        slong coeff_bits,  exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_univar_init(fx, ctx);       

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);

        exp_bits1 = n_randint(state, 3*FLINT_BITS/2) + 1;
        exp_bits2 = n_randint(state, 3*FLINT_BITS/2) + 1;
        coeff_bits = n_randint(state, 100);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {

            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);

            success = fmpz_mpoly_to_univar(fx, f, j, ctx);
            if (!success) {
                continue;
            }

            fmpz_mpoly_from_univar(g, fx, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

            if (!fmpz_mpoly_equal(f,g,ctx))
            {
                printf("FAIL\n");
                flint_printf("Check mpoly -> mpoly_univar -> mpoly\ni: %wd  j: %wd\n",i,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);  
        fmpz_mpoly_clear(g, ctx);  
        fmpz_mpoly_univar_clear(fx, ctx);       
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

