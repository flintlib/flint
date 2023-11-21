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
#include "long_extras.h"

int
arf_add_si_naive(arf_t z, const arf_t x, slong y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    int r;
    arf_init(t);
    arf_set_si(t, y);
    r = arf_add(z, x, t, prec, rnd);
    arf_clear(t);
    return r;
}

TEST_FUNCTION_START(arf_add_si, state)
{
    slong iter, iter2;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, z, v;
        slong y;
        slong prec, r1, r2;
        arf_rnd_t rnd;

        arf_init(x);
        arf_init(z);
        arf_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arf_randtest_special(x, state, 2000, 10);
            y = z_randtest(state);
            prec = 2 + n_randint(state, 2000);

            if (n_randint(state, 10) == 0 && fmpz_bits(ARF_EXPREF(x)) < 10)
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

            switch (n_randint(state, 1))
            {
            case 0:
                r1 = arf_add_si(z, x, y, prec, rnd);
                r2 = arf_add_si_naive(v, x, y, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL! (ha iter %wd, %wd)\n", iter, iter2);
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); flint_printf("%wd", y); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            default:
                r2 = arf_add_si_naive(v, x, y, prec, rnd);
                r1 = arf_add_si(x, x, y, prec, rnd);
                if (!arf_equal(x, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); flint_printf("%wd", y); flint_printf("\n\n");
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
    }

    TEST_FUNCTION_END(state);
}
