/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
arf_sosq_naive(arf_t z, const arf_t x, const arf_t y, slong prec, arf_rnd_t rnd)
{
    arf_t t;
    int inexact;
    arf_init(t);
    arf_mul(t, x, x, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul(z, y, y, ARF_PREC_EXACT, ARF_RND_DOWN);
    inexact = arf_add(z, z, t, prec, rnd);
    arf_clear(t);
    return inexact;
}

int main()
{
    slong iter, iter2;
    flint_rand_t state;

    flint_printf("sosq....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        arf_t x, y, z, v;
        slong prec, r1, r2;
        arf_rnd_t rnd;

        arf_init(x);
        arf_init(y);
        arf_init(z);
        arf_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arf_randtest_special(x, state, 2000, 100);
            arf_randtest_special(y, state, 2000, 100);
            arf_randtest_special(z, state, 2000, 100);
            arf_randtest_special(v, state, 2000, 100);

            prec = 2 + n_randint(state, 2000);

            if (n_randint(state, 10) == 0 &&
                fmpz_bits(ARF_EXPREF(x)) < 10 &&
                fmpz_bits(ARF_EXPREF(y)) < 10)
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

            switch (n_randint(state, 5))
            {
            case 0:
                r1 = arf_sosq(z, x, y, prec, rnd);
                r2 = arf_sosq_naive(v, x, y, prec, rnd);
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
                r1 = arf_sosq(z, x, x, prec, rnd);
                r2 = arf_sosq_naive(v, x, x, prec, rnd);
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
                r2 = arf_sosq_naive(v, x, y, prec, rnd);
                r1 = arf_sosq(x, x, y, prec, rnd);
                if (!arf_equal(v, x) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 2)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    flint_abort();
                }
                break;

            case 3:
                r2 = arf_sosq_naive(v, x, y, prec, rnd);
                r1 = arf_sosq(y, x, y, prec, rnd);
                if (!arf_equal(v, y) || r1 != r2)
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
                r2 = arf_sosq_naive(v, x, x, prec, rnd);
                r1 = arf_sosq(x, x, x, prec, rnd);
                if (!arf_equal(v, x) || r1 != r2)
                {
                    flint_printf("FAIL (aliasing 4)!\n");
                    flint_printf("prec = %wd, rnd = %d\n\n", prec, rnd);
                    flint_printf("x = "); arf_print(x); flint_printf("\n\n");
                    flint_printf("v = "); arf_print(v); flint_printf("\n\n");
                    flint_printf("z = "); arf_print(z); flint_printf("\n\n");
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
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
