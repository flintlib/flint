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

TEST_FUNCTION_START(acb_mul_naive, state)
{
    slong iter;

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

        acb_mul_naive(c, a, b, prec);
        acb_mul_naive(a, a, b, prec);

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

        acb_mul_naive(c, a, b, prec);
        acb_mul_naive(b, a, b, prec);

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
        acb_mul_naive(c, a, a, prec);
        acb_mul_naive(d, a, b, prec);

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
        acb_mul_naive(c, a, b, prec);
        acb_mul_naive(a, a, a, prec);

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
        acb_mul_naive(e, a, d, 2 + n_randint(state, 200));

        acb_mul_naive(d, a, b, 2 + n_randint(state, 200));
        acb_mul_naive(f, a, c, 2 + n_randint(state, 200));
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

    TEST_FUNCTION_END(state);
}
