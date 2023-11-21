/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

void
arb_fma_naive(arb_t res, const arb_t x, const arb_t y, const arb_t z, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set(t, z);
    arb_addmul(t, x, y, prec);
    arb_set(res, t);
    arb_clear(t);
}

TEST_FUNCTION_START(arb_fma, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z, res1, res2;
        slong prec;
        int aliasing;

        arb_init(x);
        arb_init(y);
        arb_init(z);
        arb_init(res1);
        arb_init(res2);

        prec = 2 + n_randint(state, 200);

        arb_randtest_special(x, state, 200, 100);
        arb_randtest_special(y, state, 200, 100);
        arb_randtest_special(z, state, 200, 100);
        arb_randtest_special(res1, state, 200, 100);
        arb_randtest_special(res2, state, 200, 100);

        if (n_randint(state, 10) == 0 &&
            fmpz_bits(ARF_EXPREF(arb_midref(x))) < 10 &&
            fmpz_bits(ARF_EXPREF(arb_midref(y))) < 10 &&
            fmpz_bits(ARF_EXPREF(arb_midref(z))) < 10)
        {
            prec = ARF_PREC_EXACT;
        }

        aliasing = n_randint(state, 7);

        switch (aliasing)
        {
            case 0:
                arb_fma(res1, x, y, z, prec);
                arb_fma_naive(res2, x, y, z, prec);
                break;
            case 1:
                arb_set(res1, z);
                arb_fma(res1, x, y, res1, prec);
                arb_fma_naive(res2, x, y, z, prec);
                break;
            case 2:
                arb_set(res1, x);
                arb_fma(res1, res1, y, z, prec);
                arb_fma_naive(res2, x, y, z, prec);
                break;
            case 3:
                arb_set(res1, y);
                arb_fma(res1, x, res1, z, prec);
                arb_fma_naive(res2, x, y, z, prec);
                break;
            case 4:
                arb_fma(res1, x, x, z, prec);
                arb_fma_naive(res2, x, x, z, prec);
                break;
            case 5:
                arb_set(res1, x);
                arb_fma(res1, res1, res1, z, prec);
                arb_fma_naive(res2, x, x, z, prec);
                break;
            default:
                arb_set(res1, x);
                arb_fma(res1, res1, res1, res1, prec);
                arb_fma_naive(res2, x, x, x, prec);
                break;
        }

        if (!arb_equal(res1, res2))
        {
            flint_printf("FAIL!\n");
            flint_printf("prec = %wd, aliasing = %d\n\n", prec, aliasing);
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 30); flint_printf("\n\n");
            flint_printf("z = "); arb_printd(z, 30); flint_printf("\n\n");
            flint_printf("res1 = "); arb_printd(res1, 30); flint_printf("\n\n");
            flint_printf("res2 = "); arb_printd(res2, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
        arb_clear(res1);
        arb_clear(res2);
    }

    TEST_FUNCTION_END(state);
}
