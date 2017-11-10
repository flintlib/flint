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
    int i, j;

    FLINT_TEST_INIT(state);

    flint_printf("univar....");
    fflush(stdout);


    /* Check mpoly -> mpoly_univar -> mpoly */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g;
        fmpz_mpoly_univar_t fx;
        ordering_t ord;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits,  exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_univar_init(fx, ctx);       

        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        exp_bits1 = n_randint(state, (FLINT_BITS - 1)/(nvars
                          + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bits2 = n_randint(state, (FLINT_BITS - 1)/(nvars
                          + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < nvars; j++)
        {

            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);

            fmpz_mpoly_to_univar(fx, f, j, ctx);
            fmpz_mpoly_from_univar(g, fx, ctx);

            if (!fmpz_mpoly_equal(f,g,ctx))
            {
                printf("FAIL\n");
                flint_printf("Check mpoly -> mpoly_univar -> mpoly\n"
                                                       "i: %wd  j: %wd\n",i,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);  
        fmpz_mpoly_clear(g, ctx);  
        fmpz_mpoly_univar_clear(fx, ctx);       
    }


    /* Check addition commutes */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h1, h2;
        fmpz_mpoly_univar_t fx, gx, h1x, h2x;
        ordering_t ord;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h1, ctx);
        fmpz_mpoly_init(h2, ctx);
        fmpz_mpoly_univar_init(fx, ctx);
        fmpz_mpoly_univar_init(gx, ctx);
        fmpz_mpoly_univar_init(h1x, ctx);
        fmpz_mpoly_univar_init(h2x, ctx);

        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);
        exp_bits1 = n_randint(state, (FLINT_BITS - 1)/(nvars
                          + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bits2 = n_randint(state, (FLINT_BITS - 1)/(nvars
                          + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < nvars; j++)
        {
            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);

            fmpz_mpoly_to_univar(fx, f, j, ctx);
            fmpz_mpoly_to_univar(gx, g, j, ctx);
            fmpz_mpoly_univar_add(h1x, fx, gx, ctx);
            fmpz_mpoly_add(h2, f, g, ctx);
            fmpz_mpoly_from_univar(h1, h1x, ctx);
            fmpz_mpoly_to_univar(h2x, h2, j, ctx);
            fmpz_mpoly_univar_add(fx, fx, gx, ctx);

            if (   !fmpz_mpoly_equal(h1, h2, ctx)
                || !fmpz_mpoly_univar_equal(h1x, h2x, ctx) 
                || !fmpz_mpoly_univar_equal(h1x, fx, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check addition commutes\ni: %wd  j: %wd\n",i,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h1, ctx);
        fmpz_mpoly_clear(h2, ctx);
        fmpz_mpoly_univar_clear(fx, ctx);       
        fmpz_mpoly_univar_clear(gx, ctx);       
        fmpz_mpoly_univar_clear(h1x, ctx);       
        fmpz_mpoly_univar_clear(h2x, ctx);       
    }

    /* Check multiplication commutes */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h1, h2;
        fmpz_mpoly_univar_t fx, gx, h1x, h2x;
        ordering_t ord;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h1, ctx);
        fmpz_mpoly_init(h2, ctx);
        fmpz_mpoly_univar_init(fx, ctx);
        fmpz_mpoly_univar_init(gx, ctx);
        fmpz_mpoly_univar_init(h1x, ctx);
        fmpz_mpoly_univar_init(h2x, ctx);

        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);
        exp_bits1 = n_randint(state, (FLINT_BITS - 3)/(nvars
                          + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bits2 = n_randint(state, (FLINT_BITS - 3)/(nvars
                          + mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < nvars; j++)
        {
            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);

            fmpz_mpoly_to_univar(fx, f, j, ctx);
            fmpz_mpoly_to_univar(gx, g, j, ctx);
            fmpz_mpoly_univar_mul(h1x, fx, gx, ctx);
            fmpz_mpoly_mul_johnson(h2, f, g, ctx);
            fmpz_mpoly_from_univar(h1, h1x, ctx);
            fmpz_mpoly_to_univar(h2x, h2, j, ctx);
            fmpz_mpoly_univar_mul(fx, fx, gx, ctx);

            if (   !fmpz_mpoly_equal(h1, h2, ctx)
                || !fmpz_mpoly_univar_equal(h1x, h2x, ctx)
                || !fmpz_mpoly_univar_equal(h1x, fx, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check multiplication commutes\n"
                                                       "i: %wd  j: %wd\n",i,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h1, ctx);
        fmpz_mpoly_clear(h2, ctx);
        fmpz_mpoly_univar_clear(fx, ctx);       
        fmpz_mpoly_univar_clear(gx, ctx);       
        fmpz_mpoly_univar_clear(h1x, ctx);       
        fmpz_mpoly_univar_clear(h2x, ctx);       
    }


    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

