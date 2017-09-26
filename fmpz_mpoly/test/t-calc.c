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
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("differentiate/integrate....");
    fflush(stdout);

    /* randomized testing doesn't catch exponent overflow in integrate */
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g;
        fmpz_t s;
        const char* vars[] = {"x","y","z","w","s","t","u","v"};

        fmpz_mpoly_ctx_init(ctx, 7, ORD_DEGLEX);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_init(s);

        result = fmpz_mpoly_set_str_pretty(f,
                           "-2*x^2+x*(-(x+y)+x^3*y^2+y^3)+z*y^126", vars, ctx);

        if (result) {
                printf("FAIL\n");
                flint_printf("set_str_pretty\n");
                flint_abort();
        }

        fmpz_mpoly_integrate(g, s, f, 1, ctx);
        fmpz_mpoly_differentiate(g, g, 1, ctx);
        fmpz_mpoly_scalar_divexact_fmpz(g, g, s, ctx);

        result = fmpz_mpoly_equal(f, f, ctx);
        if (!result) {
                printf("FAIL\n");
                flint_printf("manual integrate check\n");
                flint_abort();
        }

        fmpz_clear(s);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
    }



    /* Check d(f*g) = df*g + f*dg */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, fp, gp, hp, t1, t2;
        ordering_t ord;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong idx, coeff_bits, exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);
        fmpz_mpoly_init(hp, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, FLINT_BITS - 1 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
        exp_bits1 = n_randint(state, FLINT_BITS - 2 -
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
        exp_bits2 = n_randint(state, FLINT_BITS - 2 -
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 200);

        fmpz_mpoly_randtest(hp, state, len, exp_bound, coeff_bits, ctx);
        fmpz_mpoly_randtest(fp, state, len, exp_bound, coeff_bits, ctx);
        fmpz_mpoly_randtest(gp, state, len, exp_bound, coeff_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);

            idx = n_randint(state, nvars);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_test(h, ctx);

            fmpz_mpoly_differentiate(hp, h, idx, ctx);
            fmpz_mpoly_test(hp, ctx);
            fmpz_mpoly_differentiate(fp, f, idx, ctx);
            fmpz_mpoly_test(fp, ctx);
            fmpz_mpoly_differentiate(gp, g, idx, ctx);
            fmpz_mpoly_test(gp, ctx);

            fmpz_mpoly_mul_johnson(t1, f, gp, ctx);
            fmpz_mpoly_test(t1, ctx);
            fmpz_mpoly_mul_johnson(t2, g, fp, ctx);
            fmpz_mpoly_test(t2, ctx);
            fmpz_mpoly_add(t1, t1, t2, ctx);
            fmpz_mpoly_test(t1, ctx);

            result = fmpz_mpoly_equal(hp, t1, ctx);

            if (!result) {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                flint_abort();
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);  
        fmpz_mpoly_clear(g, ctx);  
        fmpz_mpoly_clear(h, ctx);  
        fmpz_mpoly_clear(fp, ctx);  
        fmpz_mpoly_clear(gp, ctx);  
        fmpz_mpoly_clear(hp, ctx);  
        fmpz_mpoly_clear(t1, ctx);  
        fmpz_mpoly_clear(t2, ctx);  
    }


    /* Check d(f*g) = df*g + f*dg with aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, fp, gp, t1, t2;
        ordering_t ord;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong idx, coeff_bits, exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, FLINT_BITS - 1 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
        exp_bits1 = n_randint(state, FLINT_BITS - 2 -
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
        exp_bits2 = n_randint(state, FLINT_BITS - 2 -
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 200);

        fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);
        fmpz_mpoly_randtest(fp, state, len, exp_bound, coeff_bits, ctx);
        fmpz_mpoly_randtest(gp, state, len, exp_bound, coeff_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);

            idx = n_randint(state, nvars);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_test(h, ctx);

            fmpz_mpoly_differentiate(h, h, idx, ctx);
            fmpz_mpoly_test(h, ctx);
            fmpz_mpoly_set(fp, f, ctx);
            fmpz_mpoly_differentiate(fp, fp, idx, ctx);
            fmpz_mpoly_test(fp, ctx);
            fmpz_mpoly_set(gp, g, ctx);
            fmpz_mpoly_differentiate(gp, gp, idx, ctx);
            fmpz_mpoly_test(gp, ctx);

            fmpz_mpoly_mul_johnson(t1, f, gp, ctx);
            fmpz_mpoly_test(t1, ctx);
            fmpz_mpoly_mul_johnson(t2, g, fp, ctx);
            fmpz_mpoly_test(t2, ctx);
            fmpz_mpoly_add(t1, t1, t2, ctx);
            fmpz_mpoly_test(t1, ctx);

            result = fmpz_mpoly_equal(h, t1, ctx);

            if (!result) {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg with aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                flint_abort();
           }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
       fmpz_mpoly_clear(fp, ctx);  
       fmpz_mpoly_clear(gp, ctx); 
       fmpz_mpoly_clear(t1, ctx);  
       fmpz_mpoly_clear(t2, ctx);  
    }

    /* Check d(int(f)) = f with aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, f1;
        ordering_t ord;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong idx, coeff_bits, exp_bits, exp_bits1, exp_bits2;
        fmpz_t s, s1;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_init(s);
        fmpz_init(s1);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(f1, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, FLINT_BITS - 1 - 
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
        exp_bits1 = n_randint(state, FLINT_BITS - 2 -
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;
        exp_bits2 = n_randint(state, FLINT_BITS - 2 -
                  mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars)) + 1;

        exp_bound = n_randbits(state, exp_bits);
        exp_bound1 = n_randbits(state, exp_bits1);
        exp_bound2 = n_randbits(state, exp_bits2);

        coeff_bits = n_randint(state, 200);

        fmpz_mpoly_randtest(f, state, len, exp_bound, coeff_bits, ctx);
        fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits, ctx);
        fmpz_mpoly_randtest(f1, state, len1, exp_bound1, coeff_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            idx = n_randint(state, nvars);

            fmpz_mpoly_set(g, f, ctx);
            fmpz_mpoly_test(g, ctx);
            fmpz_mpoly_integrate(f1, s1, f, idx, ctx);
            fmpz_mpoly_test(f1, ctx);
            fmpz_mpoly_integrate(f, s, f, idx, ctx);
            fmpz_mpoly_test(f, ctx);
            result = fmpz_equal(s, s1) && fmpz_mpoly_equal(f, f1, ctx);
            fmpz_mpoly_differentiate(f, f, idx, ctx);
            fmpz_mpoly_test(f, ctx);
            fmpz_mpoly_scalar_mul_fmpz(g, g, s, ctx);
            result = result && fmpz_mpoly_equal(f, g, ctx);
                                       
            if (!result) {
                printf("FAIL\n");
                flint_printf("Check d(int(f)) = f with aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                flint_abort();
            }
        }

        fmpz_clear(s);
        fmpz_clear(s1);
        fmpz_mpoly_clear(f, ctx);  
        fmpz_mpoly_clear(g, ctx);  
        fmpz_mpoly_clear(f1, ctx);  
    }


    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}

