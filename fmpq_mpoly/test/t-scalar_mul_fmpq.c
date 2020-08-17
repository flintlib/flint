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
#include "fmpq_mpoly.h"

int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("scalar_mul_fmpq....");
    fflush(stdout);

    /* Check (f*a)*b = f*(a*b) */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k;
        ordering_t ord;
        fmpq_t a, b, c;
        slong nvars, len, coeff_bits, exp_bits;

        fmpq_init(a);
        fmpq_init(b);
        fmpq_init(c);

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpq_mpoly_ctx_init(ctx, nvars, ord);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 2;
        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 10; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

            fmpq_randtest(a, state, coeff_bits + 1);
            fmpq_randtest(b, state, coeff_bits + 1);
            fmpq_mul(c, a, b);

            fmpq_mpoly_scalar_mul_fmpq(g, f, a, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);
            fmpq_mpoly_scalar_mul_fmpq(h, g, b, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_scalar_mul_fmpq(k, f, c, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);

            result = fmpq_mpoly_equal(h, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check (f*a)*b = f*(a*b)\ni=%wd j=%wd\n",i,j);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);  
        fmpq_mpoly_clear(g, ctx);  
        fmpq_mpoly_clear(h, ctx);  
        fmpq_mpoly_clear(k, ctx);  

        fmpq_clear(a);      
        fmpq_clear(b);      
        fmpq_clear(c);      
    }

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        ordering_t ord;
        fmpq_t c;
        slong nvars, len, coeff_bits, exp_bits;

        fmpq_init(c);

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 20) + 1;

        fmpq_mpoly_ctx_init(ctx, nvars, ord);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 2;
        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 10; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpq_randtest(c, state, coeff_bits + 1);

            fmpq_mpoly_set(g, f, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);
            fmpq_mpoly_scalar_mul_fmpq(h, f, c, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_scalar_mul_fmpq(g, g, c, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);
            result = fmpq_mpoly_equal(g, h, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing\ni=%wd j=%wd\n",i,j);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);

        fmpq_clear(c);      
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

