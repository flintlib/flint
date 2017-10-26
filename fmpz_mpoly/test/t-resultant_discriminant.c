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

    flint_printf("resultant_discriminant....\n");
    fflush(stdout);


    {



        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c;
        fmpz_mpoly_univar_t ax, bx, cx;
        timeit_t time;
        const char* vars[] = {"x","y","z","a","b","c","d","e","f","g","h"};

        fmpz_mpoly_ctx_init(ctx, 2, ORD_DEGREVLEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);

        fmpz_mpoly_set_str_pretty(a, "-x^4*y-2*x^4-y^3-x^2", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "2*y^2", vars, ctx);

        printf("A: "); fmpz_mpoly_print_pretty(a, vars, ctx); printf("\n");
        printf("B: "); fmpz_mpoly_print_pretty(b, vars, ctx); printf("\n");
        fmpz_mpoly_resultant(c, a, b, 1, ctx);
        printf("C: "); fmpz_mpoly_print_pretty(c, vars, ctx); printf("\n");

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
    }



    /* Check res(a*b,c) = res(a,c)*res(b,c) */
    for (i = 0; i < 0*30 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c, ab, ra, rb, rab, p;
        fmpz_mpoly_univar_t fx;
        ordering_t ord;
        slong nvars, len1, len2, len3, exp_bound, exp_bound1, exp_bound2, exp_bound3;
        slong coeff_bits, exp_bits1, exp_bits2, exp_bits3;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 2) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);       
        fmpz_mpoly_init(ab, ctx);       
        fmpz_mpoly_init(ra, ctx);
        fmpz_mpoly_init(rb, ctx);
        fmpz_mpoly_init(rab, ctx);
        fmpz_mpoly_init(p, ctx);

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);
        len3 = n_randint(state, 20);
        exp_bound1 = n_randint(state, 8) + 1;
        exp_bound2 = n_randint(state, 8) + 1;
        exp_bound3 = n_randint(state, 8) + 1;
        coeff_bits = n_randint(state, 3);
        fmpz_mpoly_randtest(a, state, len1, exp_bound1, coeff_bits, ctx);
        fmpz_mpoly_randtest(b, state, len2, exp_bound2, coeff_bits, ctx);
        fmpz_mpoly_randtest(c, state, len3, exp_bound3, coeff_bits, ctx);

        for (j = 0; j < nvars; j++)
        {
            fmpz_mpoly_mul_johnson(ab, a, b, ctx);
            fmpz_mpoly_resultant(ra, a, c, j, ctx);
            fmpz_mpoly_resultant(rb, b, c, j, ctx);
            fmpz_mpoly_resultant(rab, ab, c, j, ctx);
            fmpz_mpoly_mul_johnson(p, ra, rb, ctx);

            if (!fmpz_mpoly_equal(p,rab,ctx))
            {
                printf("FAIL\n");
                flint_printf(" Check res(a*b,c) = res(a,c)*res(b,c) \ni: %wd  j: %wd\n",i,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);       
        fmpz_mpoly_clear(ab, ctx);       
        fmpz_mpoly_clear(ra, ctx);
        fmpz_mpoly_clear(rb, ctx);
        fmpz_mpoly_clear(rab, ctx);
        fmpz_mpoly_clear(p, ctx);
    }


    /* Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2 */
    for (i = 0; i < 0*30 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, ab, r, da, db, dab, p;
        ordering_t ord;
        slong nvars, len1, len2, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 2) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(ab, ctx);
        fmpz_mpoly_init(da, ctx);
        fmpz_mpoly_init(db, ctx);
        fmpz_mpoly_init(dab, ctx);
        fmpz_mpoly_init(r, ctx);
        fmpz_mpoly_init(p, ctx);

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20);
        exp_bound1 = n_randint(state, 8) + 1;
        exp_bound2 = n_randint(state, 8) + 1;
        coeff_bits = n_randint(state, 3);
        fmpz_mpoly_randtest(a, state, len1, exp_bound1, coeff_bits, ctx);
        fmpz_mpoly_randtest(b, state, len2, exp_bound2, coeff_bits, ctx);

        for (j = 0; j < nvars; j++)
        {
            if (fmpz_mpoly_degree(a, j, ctx) < 1)
                continue;
            if (fmpz_mpoly_degree(b, j, ctx) < 1)
                continue;

            fmpz_mpoly_mul_johnson(ab, a, b, ctx);
            fmpz_mpoly_resultant(r, a, b, j, ctx);
            fmpz_mpoly_discriminant(da, a, j, ctx);
            fmpz_mpoly_discriminant(db, b, j, ctx);
            fmpz_mpoly_discriminant(dab, ab, j, ctx);
            fmpz_mpoly_mul_johnson(p, da, db, ctx);
            fmpz_mpoly_mul_johnson(p, p, r, ctx);
            fmpz_mpoly_mul_johnson(p, p, r, ctx);

            if (!fmpz_mpoly_equal(dab, p, ctx))
            {
                printf("FAIL\n");
                flint_printf(" Check disc(a*b) = disc(a)*disc(b)*res(a,b)^2 \ni: %wd  j: %wd\n",i,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(ab, ctx);       
        fmpz_mpoly_clear(da, ctx);
        fmpz_mpoly_clear(db, ctx);
        fmpz_mpoly_clear(dab, ctx);
        fmpz_mpoly_clear(r, ctx);
        fmpz_mpoly_clear(p, ctx);
    }



    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

