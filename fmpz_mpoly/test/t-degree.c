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
#include "profiler.h"

int
main(void)
{
    int i, j;

    FLINT_TEST_INIT(state);

    flint_printf("degree....");
    fflush(stdout);


    /* Check degree does not go up under addition */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        ordering_t ord;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        exp_bits1 = n_randint(state, (FLINT_BITS - 1)/(nvars
                          + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bits2 = n_randint(state, (FLINT_BITS - 1)/(nvars
                          + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 4);

        for (j = 0; j < nvars; j++)
        {
            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);
            fmpz_mpoly_add(h, f, g, ctx);

            if (fmpz_mpoly_degree(h, j, ctx) > FLINT_MAX(
                                                  fmpz_mpoly_degree(f, j, ctx),
                                                  fmpz_mpoly_degree(g, j, ctx)
                                               )
               )
            {
                printf("FAIL\n");
                flint_printf("Check degree does not go up under addition\n"
                                                     "i: %wd  j: %wd\n", i, j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx); 
    }


    /* Check degree adds under multiplication */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        ordering_t ord;
        slong df, dg, dh;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        exp_bits1 = n_randint(state, (FLINT_BITS - 1)/(nvars + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bits2 = n_randint(state, (FLINT_BITS - 1)/(nvars + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 4);

        for (j = 0; j < nvars; j++)
        {
            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);
            fmpz_mpoly_mul_johnson(h, f, g, ctx);

            df = fmpz_mpoly_degree(f, j, ctx);
            dg = fmpz_mpoly_degree(g, j, ctx);
            dh = fmpz_mpoly_degree(h, j, ctx);



            if (dh != ((df + dg) | FLINT_SIGN_EXT(df) | FLINT_SIGN_EXT(dg)))
            {
                printf("FAIL\n");
                flint_printf("Check degree adds under multiplication\n"
                                                     "i: %wd  j: %wd\n", i, j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx); 
    }


    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

