/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arf.h"

TEST_FUNCTION_START(arf_complex_sqr, state)
{
    slong iter, iter2;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arf_t e1, f1, e2, f2, a, b;
        slong prec, r1, r2;
        arf_rnd_t rnd;

        arf_init(a);
        arf_init(b);
        arf_init(e1);
        arf_init(f1);
        arf_init(e2);
        arf_init(f2);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            if (n_randint(state, 100) == 0)
            {
                arf_randtest_special(a, state, 12000, 3);
                arf_randtest_special(b, state, 12000, 3);
                prec = 2 + n_randint(state, 12000);
            }
            else
            {
                arf_randtest_special(a, state, 3000, 100);
                arf_randtest_special(b, state, 3000, 100);
                prec = 2 + n_randint(state, 3000);
            }

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
                r1 = arf_complex_sqr(e1, f1, a, b, prec, rnd);
                r2 = arf_complex_mul_fallback(e2, f2, a, b, a, b, prec, rnd);
                if (!arf_equal(e1, e2) || !arf_equal(f1, f2) || r1 != r2)
                {
                    flint_printf("FAIL!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("a = "); arf_print(a); flint_printf("\n\n");
                    flint_printf("b = "); arf_print(b); flint_printf("\n\n");
                    flint_printf("e1 = "); arf_print(e1); flint_printf("\n\n");
                    flint_printf("f1 = "); arf_print(f1); flint_printf("\n\n");
                    flint_printf("e2 = "); arf_print(e2); flint_printf("\n\n");
                    flint_printf("f2 = "); arf_print(f2); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            default:
                r1 = arf_complex_mul_fallback(e1, f1, a, b, a, b, prec, rnd);
                r2 = arf_complex_sqr(a, b, a, b, prec, rnd);
                if (!arf_equal(e1, a) || !arf_equal(f1, b) || r1 != r2)
                {
                    flint_printf("FAIL! (aliasing)\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("a = "); arf_print(a); flint_printf("\n\n");
                    flint_printf("b = "); arf_print(b); flint_printf("\n\n");
                    flint_printf("e1 = "); arf_print(e1); flint_printf("\n\n");
                    flint_printf("f1 = "); arf_print(f1); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;
            }
        }

        arf_clear(a);
        arf_clear(b);
        arf_clear(e1);
        arf_clear(f1);
        arf_clear(e2);
        arf_clear(f2);
    }

    TEST_FUNCTION_END(state);
}
