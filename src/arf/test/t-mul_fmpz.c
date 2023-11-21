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
arf_mul_fmpz_naive(arf_t z, const arf_t x, const fmpz_t y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    int r;
    arf_init(t);
    arf_set_fmpz(t, y);
    r = arf_mul(z, x, t, prec, rnd);
    arf_clear(t);
    return r;
}

TEST_FUNCTION_START(arf_mul_fmpz, state)
{
    slong iter, iter2;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, z, v;
        fmpz_t y;
        slong prec, r1, r2;
        arf_rnd_t rnd;

        arf_init(x);
        arf_init(z);
        arf_init(v);
        fmpz_init(y);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arf_randtest_special(x, state, 2000, 10);
            fmpz_randtest(y, state, 2000);
            prec = 2 + n_randint(state, 2000);

            if (n_randint(state, 10))
                prec = ARF_PREC_EXACT;

            switch (n_randint(state, 5))
            {
                case 0:  rnd = ARF_RND_DOWN; break;
                case 1:  rnd = ARF_RND_UP; break;
                case 2:  rnd = ARF_RND_FLOOR; break;
                case 3:  rnd = ARF_RND_CEIL; break;
                default: rnd = ARF_RND_NEAR; break;
            }

            switch (n_randint(state, 2))
            {
            case 0:
                r1 = arf_mul_fmpz(z, x, y, prec, rnd);
                r2 = arf_mul_fmpz_naive(v, x, y, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            default:
                r2 = arf_mul_fmpz_naive(v, x, y, prec, rnd);
                r1 = arf_mul_fmpz(x, x, y, prec, rnd);
                if (!arf_equal(x, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;
            }
        }

        arf_clear(x);
        arf_clear(z);
        arf_clear(v);
        fmpz_clear(y);
    }

    TEST_FUNCTION_END(state);
}
