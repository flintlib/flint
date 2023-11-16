/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

/* these functions are not exposed to the public for now,
   but it still makes sense to test them explicitly */
void arb_exp_taylor_sum_rs_generic(arb_t s, const arb_t x, slong N, slong prec);

TEST_FUNCTION_START(arb_exp_arf_rs_generic, state)
{
    slong iter;

    /* test the rs algorithm explicitly */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, s1, s2;
        slong prec;
        int minus1;

        arb_init(x);
        arb_init(s1);
        arb_init(s2);

        prec = 2 + n_randint(state, 2000);
        minus1 = n_randint(state, 2);
        arb_randtest(x, state, 1 + n_randint(state, 2000), 3);
        mag_zero(arb_radref(x));

        if (n_randint(state, 2))
            arb_mul_2exp_si(x, x, -n_randint(state, 2 * prec));

        if (minus1)
            arb_expm1(s1, x, prec);
        else
            arb_exp(s1, x, prec);

        switch (n_randint(state, 2))
        {
            case 0:
                arb_exp_arf_rs_generic(s2, arb_midref(x), prec, minus1);
                break;
            case 1:
                arb_set(s2, x);
                arb_exp_arf_rs_generic(s2, arb_midref(s2), prec, minus1);
                break;
        }

        if (!arb_overlaps(s1, s2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("prec = %wd, minus1 = %d\n\n", prec, minus1);
            flint_printf("x = "); arb_printn(x, 500, 0); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 500, 0); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 500, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_rel_accuracy_bits(s2) < prec - 2)
        {
            flint_printf("FAIL: poor accuracy\n\n");
            flint_printf("prec = %wd,  acc1 = %wd,  minus1 = %d\n\n",
                prec, arb_rel_accuracy_bits(s2), minus1);
            flint_printf("x = "); arb_printn(x, 500, 0); flint_printf("\n\n");
            flint_printf("s1 = "); arb_printn(s1, 500, 0); flint_printf("\n\n");
            flint_printf("s2 = "); arb_printn(s2, 500, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(s1);
        arb_clear(s2);
    }

    /* test the series evaluation code directly */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z;
        slong prec;
        slong N;

        arb_init(x);
        arb_init(y);
        arb_init(z);

        prec = 2 + n_randint(state, 2000);
        N = n_randint(state, 100);
        arb_randtest(x, state, 1 + n_randint(state, 2000), 1);
        mag_zero(arb_radref(x));

        if (n_randint(state, 2))
            arb_mul_2exp_si(x, x, -n_randint(state, 2 * prec));

        arb_exp(y, x, prec);
        arb_exp_taylor_sum_rs_generic(z, x, N, prec);

        if (!arb_overlaps(z, y))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("prec = %wd, N = %wd\n\n", prec, N);
            flint_printf("x = "); arb_printn(x, 500, 0); flint_printf("\n\n");
            flint_printf("y = "); arb_printn(y, 500, 0); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 500, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
    }

    TEST_FUNCTION_END(state);
}
