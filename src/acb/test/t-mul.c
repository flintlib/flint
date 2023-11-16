/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_mul, state)
{
    slong iter, iter2;

    /* test aliasing of c and a */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c;
        slong prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_randtest(a, state, 1 + n_randint(state, 200), 10);
        acb_randtest(b, state, 1 + n_randint(state, 200), 10);
        acb_randtest(c, state, 1 + n_randint(state, 200), 10);

        prec = 2 + n_randint(state, 200);

        acb_mul(c, a, b, prec);
        acb_mul(a, a, b, prec);

        if (!acb_equal(a, c))
        {
            flint_printf("FAIL: aliasing c, a\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }

    /* test aliasing of c and b */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c;
        slong prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_randtest(a, state, 1 + n_randint(state, 200), 10);
        acb_randtest(b, state, 1 + n_randint(state, 200), 10);
        acb_randtest(c, state, 1 + n_randint(state, 200), 10);

        prec = 2 + n_randint(state, 200);

        acb_mul(c, a, b, prec);
        acb_mul(b, a, b, prec);

        if (!acb_equal(b, c))
        {
            flint_printf("FAIL: aliasing b, a\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }

    /* test aliasing a, a */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, d;
        slong prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);

        acb_randtest(a, state, 1 + n_randint(state, 200), 10);
        acb_randtest(b, state, 1 + n_randint(state, 200), 10);
        acb_randtest(c, state, 1 + n_randint(state, 200), 10);

        prec = 2 + n_randint(state, 200);

        acb_set(b, a);
        acb_mul(c, a, a, prec);
        acb_mul(d, a, b, prec);

        if (!acb_overlaps(c, d))
        {
            flint_printf("FAIL: aliasing a, a\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("d = "); acb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
    }

    /* test aliasing a, a, a */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c;
        slong prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_randtest(a, state, 1 + n_randint(state, 200), 10);
        acb_randtest(b, state, 1 + n_randint(state, 200), 10);
        acb_randtest(c, state, 1 + n_randint(state, 200), 10);

        prec = 2 + n_randint(state, 200);

        acb_set(b, a);
        acb_mul(c, a, b, prec);
        acb_mul(a, a, a, prec);

        if (!acb_overlaps(a, c))
        {
            flint_printf("FAIL: aliasing a, a, a\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }

    /* test a*(b+c) = a*b + a*c */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, d, e, f;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(d);
        acb_init(e);
        acb_init(f);

        acb_randtest(a, state, 1 + n_randint(state, 200), 10);
        acb_randtest(b, state, 1 + n_randint(state, 200), 10);
        acb_randtest(c, state, 1 + n_randint(state, 200), 10);

        acb_add(d, b, c, 2 + n_randint(state, 200));
        acb_mul(e, a, d, 2 + n_randint(state, 200));

        acb_mul(d, a, b, 2 + n_randint(state, 200));
        acb_mul(f, a, c, 2 + n_randint(state, 200));
        acb_add(f, d, f, 2 + n_randint(state, 200));

        if (!acb_overlaps(e, f))
        {
            flint_printf("FAIL: a*(b+c) = a*b + a*c\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_printf("e = "); acb_print(e); flint_printf("\n\n");
            flint_printf("f = "); acb_print(f); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(d);
        acb_clear(e);
        acb_clear(f);
    }

    /* compare with mul_naive */
    /* main test */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, y, z, v;
        slong prec;

        acb_init(x);
        acb_init(y);
        acb_init(z);
        acb_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            acb_randtest_special(x, state, n_randint(state,2) ? 2000 : 200, 200);
            acb_randtest_special(y, state, n_randint(state,2) ? 2000 : 200, 200);

            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 5))
            {
            case 0:
                acb_mul(z, x, y, prec);
                acb_mul_naive(v, x, y, prec);

                if (!acb_overlaps(z, v))
                {
                    flint_printf("FAIL!\n");
                    flint_printf("x = "); acb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); acb_print(y); flint_printf("\n\n");
                    flint_printf("z = "); acb_print(z); flint_printf("\n\n");
                    flint_printf("v = "); acb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 1:
                acb_set(y, x);
                acb_mul(z, x, y, prec);
                acb_mul(v, x, x, prec);

                if (!acb_overlaps(z, v))
                {
                    flint_printf("FAIL (aliasing 1)!\n");
                    flint_printf("x = "); acb_print(x); flint_printf("\n\n");
                    flint_printf("z = "); acb_print(z); flint_printf("\n\n");
                    flint_printf("v = "); acb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 2:
                acb_mul(v, x, x, prec);
                acb_mul(x, x, x, prec);

                if (!acb_equal(v, x))
                {
                    flint_printf("FAIL (aliasing 2)!\n");
                    flint_printf("x = "); acb_print(x); flint_printf("\n\n");
                    flint_printf("z = "); acb_print(z); flint_printf("\n\n");
                    flint_printf("v = "); acb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            case 3:
                acb_mul(v, x, y, prec);
                acb_mul(x, x, y, prec);

                if (!acb_equal(v, x))
                {
                    flint_printf("FAIL (aliasing 3)!\n");
                    flint_printf("x = "); acb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); acb_print(y); flint_printf("\n\n");
                    flint_printf("v = "); acb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;

            default:
                acb_mul(v, x, y, prec);
                acb_mul(x, y, x, prec);

                if (!acb_overlaps(v, x))
                {
                    flint_printf("FAIL (aliasing 4)!\n");
                    flint_printf("x = "); acb_print(x); flint_printf("\n\n");
                    flint_printf("y = "); acb_print(y); flint_printf("\n\n");
                    flint_printf("v = "); acb_print(v); flint_printf("\n\n");
                    flint_abort();
                }
                break;
            }
        }

        acb_clear(x);
        acb_clear(y);
        acb_clear(z);
        acb_clear(v);
    }

    TEST_FUNCTION_END(state);
}
