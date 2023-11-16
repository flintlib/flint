/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

int
arf_fma_naive(arf_t res, const arf_t x, const arf_t y, const arf_t z, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    int inexact;

    arf_init(t);
    arf_mul(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);

    inexact = arf_add(res, z, t, prec, rnd);

    arf_clear(t);

    return inexact;
}

TEST_FUNCTION_START(arf_fma, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z, res1, res2;
        slong prec, r1, r2;
        arf_rnd_t rnd;
        int aliasing;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(res1);
        arf_init(res2);

        prec = 2 + n_randint(state, 200);

        arf_randtest_special(x, state, 200, 100);
        arf_randtest_special(y, state, 200, 100);
        arf_randtest_special(z, state, 200, 100);
        arf_randtest_special(res1, state, 200, 100);
        arf_randtest_special(res2, state, 200, 100);

        if (n_randint(state, 10) == 0 &&
            fmpz_bits(ARF_EXPREF(x)) < 10 &&
            fmpz_bits(ARF_EXPREF(y)) < 10 &&
            fmpz_bits(ARF_EXPREF(z)) < 10)
        {
            prec = ARF_PREC_EXACT;
        }

        switch (n_randint(state, 5))
        {
            case 0:  rnd = ARF_RND_DOWN; break;
            case 1:  rnd = ARF_RND_UP; break;
            case 2:  rnd = ARF_RND_FLOOR; break;
            case 3:  rnd = ARF_RND_CEIL; break;
            default: rnd = ARF_RND_NEAR; break;
        }

        aliasing = n_randint(state, 7);

        switch (aliasing)
        {
            case 0:
                r1 = arf_fma(res1, x, y, z, prec, rnd);
                r2 = arf_fma_naive(res2, x, y, z, prec, rnd);
                break;
            case 1:
                arf_set(res1, z);
                r1 = arf_fma(res1, x, y, res1, prec, rnd);
                r2 = arf_fma_naive(res2, x, y, z, prec, rnd);
                break;
            case 2:
                arf_set(res1, x);
                r1 = arf_fma(res1, res1, y, z, prec, rnd);
                r2 = arf_fma_naive(res2, x, y, z, prec, rnd);
                break;
            case 3:
                arf_set(res1, y);
                r1 = arf_fma(res1, x, res1, z, prec, rnd);
                r2 = arf_fma_naive(res2, x, y, z, prec, rnd);
                break;
            case 4:
                r1 = arf_fma(res1, x, x, z, prec, rnd);
                r2 = arf_fma_naive(res2, x, x, z, prec, rnd);
                break;
            case 5:
                arf_set(res1, x);
                r1 = arf_fma(res1, res1, res1, z, prec, rnd);
                r2 = arf_fma_naive(res2, x, x, z, prec, rnd);
                break;
            default:
                arf_set(res1, x);
                r1 = arf_fma(res1, res1, res1, res1, prec, rnd);
                r2 = arf_fma_naive(res2, x, x, x, prec, rnd);
                break;
        }

        if (!arf_equal(res1, res2) || r1 != r2)
        {
            flint_printf("FAIL!\n");
            flint_printf("prec = %wd, rnd = %d, aliasing = %d\n\n", prec, rnd, aliasing);
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); arf_print(y); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            flint_printf("res1 = "); arf_print(res1); flint_printf("\n\n");
            flint_printf("res2 = "); arf_print(res2); flint_printf("\n\n");
            flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
            flint_abort();
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        arf_clear(res1);
        arf_clear(res2);
    }

    TEST_FUNCTION_END(state);
}
