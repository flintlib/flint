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
    int i, j, result;
    timeit_t time;
    long t0, t1, t2;

    FLINT_TEST_INIT(state);


    flint_printf("quasidivrem_heap....");
    fflush(stdout);

    /* Check f*g/g = f */
    t0 = t1 = t2 = 0;
    printf("\n");
    for (i = 0; i < 0 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k, r0, q0, r1, q1, r2, q2;
        fmpz_t s1, s2;
        ordering_t ord;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 10) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(q0, ctx);
        fmpz_mpoly_init(r0, ctx);
        fmpz_mpoly_init(q1, ctx);
        fmpz_mpoly_init(r1, ctx);
        fmpz_mpoly_init(q2, ctx);
        fmpz_mpoly_init(r2, ctx);
        fmpz_init(s1);
        fmpz_init(s2);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100) + 100;
        len2 = n_randint(state, 100) + 100;

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

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest(f, state, len1, exp_bound1, coeff_bits, ctx);
            do {
                fmpz_mpoly_randtest(g, state, len2, exp_bound2, coeff_bits + 1, ctx);
            } while (g->length == 0);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);

if (n_randint(state, 2)) {
            timeit_start(time);
            fmpz_mpoly_quasidivrem_heapV1(s1, q1, r1, h, g, ctx);
            timeit_stop(time);
            t1 += time->wall;
            timeit_start(time);
            fmpz_mpoly_quasidivrem_heapV2(s2, q2, r2, h, g, ctx);
            timeit_stop(time);
            t2 += time->wall;
} else {
            timeit_start(time);
            fmpz_mpoly_quasidivrem_heapV2(s2, q2, r2, h, g, ctx);
            timeit_stop(time);
            t2 += time->wall;
            timeit_start(time);
            fmpz_mpoly_quasidivrem_heapV1(s1, q1, r1, h, g, ctx);
            timeit_stop(time);
            t1 += time->wall;
}

            timeit_start(time);
            fmpz_mpoly_divrem_monagan_pearce(q0, r0, h, g, ctx);
            timeit_stop(time);
            t0 += time->wall;
            fmpz_mpoly_test(q0, ctx);
            fmpz_mpoly_test(r0, ctx);
            fmpz_mpoly_remainder_test(r1, g, ctx);

            result = fmpz_equal_ui(s1, WORD(1)) && fmpz_equal_ui(s2, WORD(1))
                                               && fmpz_mpoly_equal(q1, q2, ctx)
                                               && fmpz_mpoly_equal(q0, q1, ctx)
                                               && fmpz_mpoly_equal(q0, f, ctx);

            if (!result)
            {
                printf("FAIL1\n");
                flint_printf("i=%wd j=%wd\n",i,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(q0, ctx);
        fmpz_mpoly_clear(r0, ctx);  
        fmpz_mpoly_clear(q1, ctx);
        fmpz_mpoly_clear(r1, ctx);  
        fmpz_mpoly_clear(q2, ctx);  
        fmpz_mpoly_clear(r2, ctx);  
        fmpz_clear(s1);
        fmpz_clear(s2);

    }
    flint_printf("exact       divrem: %wd\n",t0);
    flint_printf("exact quasidivrem1: %wd\n",t1);
    flint_printf("exact quasidivrem2: %wd\n",t2);

    

    /* Check s*f = g*q + r for random polys */
    t1 = t2 = 0;
    printf("\n");

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
/*
flint_printf("i=%wd\n",i);
*/
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, q0, q1, q2, r0, r1, r2, k1, k2, l1, l2;
       fmpz_t s0, s1, s2;
       ordering_t ord;
       slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;
       slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

        do {
            ord = mpoly_ordering_randtest(state);
        } while (ord==ORD_LEX);
	   
       nvars = n_randint(state, 5) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_init(s0);
       fmpz_init(s1);
       fmpz_init(s2);
       fmpz_mpoly_init(q0, ctx);
       fmpz_mpoly_init(r0, ctx);
       fmpz_mpoly_init(q1, ctx);
       fmpz_mpoly_init(r1, ctx);
       fmpz_mpoly_init(k1, ctx);
       fmpz_mpoly_init(l1, ctx);
       fmpz_mpoly_init(q2, ctx);
       fmpz_mpoly_init(r2, ctx);
       fmpz_mpoly_init(k2, ctx);
       fmpz_mpoly_init(l2, ctx);

       len = n_randint(state, 10);
       len1 = n_randint(state, 80);
       len2 = n_randint(state, 40);

       exp_bits = n_randint(state, 10/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits1 = n_randint(state, 10/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;
       exp_bits2 = n_randint(state, 10/(nvars + 
                            mpoly_ordering_isdeg(ord) + (nvars == 1)) + 1) + 1;

       exp_bound = n_randbits(state, exp_bits);
       exp_bound1 = n_randbits(state, exp_bits1);
       exp_bound2 = n_randbits(state, exp_bits2);

       coeff_bits = n_randint(state, 50);

       for (j = 0; j < 4; j++)
       {
            fmpz_mpoly_randtest(f, state, len1+20, exp_bound1+6, coeff_bits, ctx);
            do {
                fmpz_mpoly_randtest(g, state, len2+10, exp_bound2+2, coeff_bits+10, ctx);
            } while (g->length < 2);


/*
flint_printf("numerator length: %wd  denominator length: %wd  ",f->length,g->length);
printf("leading coeff: "); fmpz_print(g->coeffs); printf("\n");
*/


            timeit_start(time);
            fmpz_mpoly_quasidivrem_heap(s0, q0, r0, f, g, ctx);
            timeit_stop(time);
            t0 += time->wall;
/*
printf("**********\n");
printf("f: "); fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n");
printf("g: "); fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n");
printf("s0: "); fmpz_print(s0); printf("\n");
printf("q0: "); fmpz_mpoly_print_pretty(q0, NULL, ctx); printf("\n");
printf("r0: "); fmpz_mpoly_print_pretty(r0, NULL, ctx); printf("\n");
*/


if (n_randint(state, 2)) {
            timeit_start(time);
            fmpz_mpoly_quasidivrem_heapV1(s1, q1, r1, f, g, ctx);
            timeit_stop(time);
            t1 += time->wall;
            timeit_start(time);
            fmpz_mpoly_quasidivrem_heapV2(s2, q2, r2, f, g, ctx);
            timeit_stop(time);
            t2 += time->wall;
} else {
            timeit_start(time);
            fmpz_mpoly_quasidivrem_heapV2(s2, q2, r2, f, g, ctx);
            timeit_stop(time);
            t2 += time->wall;
            timeit_start(time);
            fmpz_mpoly_quasidivrem_heapV1(s1, q1, r1, f, g, ctx);
            timeit_stop(time);
            t1 += time->wall;
}
/*
printf("**********\n");
printf("f: "); fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n");
printf("g: "); fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n");
printf("s1: "); fmpz_print(s1); printf("\n");
printf("q1: "); fmpz_mpoly_print_pretty(q1, NULL, ctx); printf("\n");
printf("r1: "); fmpz_mpoly_print_pretty(r1, NULL, ctx); printf("\n");

printf("**********\n");
printf("f: "); fmpz_mpoly_print_pretty(f, NULL, ctx); printf("\n");
printf("g: "); fmpz_mpoly_print_pretty(g, NULL, ctx); printf("\n");
printf("s2: "); fmpz_print(s2); printf("\n");
printf("q2: "); fmpz_mpoly_print_pretty(q2, NULL, ctx); printf("\n");
printf("r2: "); fmpz_mpoly_print_pretty(r2, NULL, ctx); printf("\n");
*/



            fmpz_mpoly_test(q2, ctx);
            fmpz_mpoly_test(r2, ctx);
            fmpz_mpoly_remainder_strongtest(r2, g, ctx);

            fmpz_mpoly_test(q1, ctx);
            fmpz_mpoly_test(r1, ctx);
            fmpz_mpoly_remainder_strongtest(r1, g, ctx);


            fmpz_mpoly_mul_johnson(k1, q1, g, ctx);
            fmpz_mpoly_add(k1, k1, r1, ctx);
            fmpz_mpoly_scalar_mul_fmpz(l1, f, s1, ctx);

            fmpz_mpoly_mul_johnson(k2, q2, g, ctx);
            fmpz_mpoly_add(k2, k2, r2, ctx);
            fmpz_mpoly_scalar_mul_fmpz(l2, f, s2, ctx);

/*
flint_printf("quotient length: %wd  remainder length: %wd\n",q1->length,r2->length);
*/
            result = fmpz_equal(s0,s1) && fmpz_mpoly_equal(q0, q1, ctx)
                                       && fmpz_mpoly_equal(r0, r1, ctx)
                                       && fmpz_mpoly_equal(l1, k1, ctx)
                                       && fmpz_mpoly_equal(l2, k2, ctx);

            if (!result)
            {
                printf("FAIL2\n");
                flint_printf("i=%wd j=%wd\n",i,j);
                flint_abort();
            }

        }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_clear(s0);
       fmpz_clear(s1);
       fmpz_clear(s2);
       fmpz_mpoly_clear(q0, ctx);
       fmpz_mpoly_clear(r0, ctx);
       fmpz_mpoly_clear(q1, ctx);
       fmpz_mpoly_clear(r1, ctx);
       fmpz_mpoly_clear(k1, ctx);
       fmpz_mpoly_clear(l1, ctx);
       fmpz_mpoly_clear(q2, ctx);
       fmpz_mpoly_clear(r2, ctx);
       fmpz_mpoly_clear(k2, ctx);
       fmpz_mpoly_clear(l2, ctx);
    }

    flint_printf("random quasidivrem0: %wd\n",t1);
    flint_printf("random quasidivrem1: %wd\n",t1);
    flint_printf("random quasidivrem2: %wd\n",t2);


    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

