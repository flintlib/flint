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
arf_sub_naive(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
{
    int r;
    arf_t t;
    arf_init(t);
    arf_neg(t, y);
    r = arf_add(z, x, t, prec, rnd);
    arf_clear(t);
    return r;
}

TEST_FUNCTION_START(arf_sub, state)
{
    slong iter, iter2;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t x, y, z, v;
        slong prec, r1, r2;
        arf_rnd_t rnd;
        fmpz_t t;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(v);
        fmpz_init(t);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arf_randtest_special(x, state, 2000, 100);
            arf_randtest_special(y, state, 2000, 100);
            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 5))
            {
                case 0:  rnd = ARF_RND_DOWN; break;
                case 1:  rnd = ARF_RND_UP; break;
                case 2:  rnd = ARF_RND_FLOOR; break;
                case 3:  rnd = ARF_RND_CEIL; break;
                default: rnd = ARF_RND_NEAR; break;
            }

            rnd = ARF_RND_DOWN;

            if (arf_is_normal(x) && arf_is_normal(y))
            {
                fmpz_sub(t, ARF_EXPREF(x), ARF_EXPREF(y));

                /* if not too far apart, sometimes test exact addition */
                if (fmpz_bits(t) < 10)
                {
                    if (n_randint(state, 10) == 0)
                        prec = ARF_PREC_EXACT;
                }
                else if (rnd == ARF_RND_NEAR)
                { /* large shift not supported in add_naive */
                    rnd = ARF_RND_DOWN;
                }
            }

            switch (n_randint(state, 5))
            {
            case 0:
                r1 = arf_sub(z, x, y, prec, rnd);
                r2 = arf_sub_naive(v, x, y, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            case 1:
                r1 = arf_sub(z, x, x, prec, rnd);
                r2 = arf_sub_naive(v, x, x, prec, rnd);
                if (!arf_equal(z, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 1)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            case 2:
                r2 = arf_sub_naive(v, x, x, prec, rnd);
                r1 = arf_sub(x, x, x, prec, rnd);
                if (!arf_equal(v, x) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 2)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            case 3:
                r2 = arf_sub_naive(v, x, y, prec, rnd);
                r1 = arf_sub(x, x, y, prec, rnd);
                if (!arf_equal(x, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 3)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            default:
                r2 = arf_sub_naive(v, y, x, prec, rnd);
                r1 = arf_sub(x, y, x, prec, rnd);
                if (!arf_equal(x, v) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 4)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("y = "); arf_print(y); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;
            }
        }

        arf_clear(x);
        arf_clear(y);
        arf_clear(z);
        arf_clear(v);
        fmpz_clear(t);
    }

    TEST_FUNCTION_END(state);
}
