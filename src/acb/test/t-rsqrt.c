/*
    Copyright (C) 2013, 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"

TEST_FUNCTION_START(acb_rsqrt, state)
{
    slong iter;

    /* check union */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, x, ra, rb, rc, rx;

        acb_init(a);
        acb_init(b);
        acb_init(c);
        acb_init(x);
        acb_init(ra);
        acb_init(rb);
        acb_init(rc);
        acb_init(rx);

        acb_randtest_precise(a, state, 1 + n_randint(state, 500), 1 + n_randint(state, 80));
        acb_randtest_precise(b, state, 1 + n_randint(state, 500), 1 + n_randint(state, 80));
        acb_randtest_precise(c, state, 1 + n_randint(state, 500), 1 + n_randint(state, 80));
        acb_randtest(ra, state, 1 + n_randint(state, 500), 1 + n_randint(state, 80));
        acb_randtest(rb, state, 1 + n_randint(state, 500), 1 + n_randint(state, 80));
        acb_randtest(rc, state, 1 + n_randint(state, 500), 1 + n_randint(state, 80));
        acb_randtest(rx, state, 1 + n_randint(state, 500), 1 + n_randint(state, 80));

        acb_union(x, a, b, 2 + n_randint(state, 500));
        acb_union(x, x, c, 2 + n_randint(state, 500));
        acb_rsqrt(rx, x, 2 + n_randint(state, 500));

        acb_rsqrt(ra, a, 2 + n_randint(state, 500));
        acb_rsqrt(rb, b, 2 + n_randint(state, 500));
        acb_rsqrt(rc, c, 2 + n_randint(state, 500));

        if (!acb_overlaps(rx, ra))
        {
            flint_printf("FAIL: overlap a\n\n");
            flint_printf("a = "); acb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); acb_printn(b, 50, 0); flint_printf("\n\n");
            flint_printf("c = "); acb_printn(c, 50, 0); flint_printf("\n\n");
            flint_printf("x = "); acb_printn(x, 50, 0); flint_printf("\n\n");
            flint_printf("ra = "); acb_printn(ra, 50, 0); flint_printf("\n\n");
            flint_printf("rx = "); acb_printn(rx, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (!acb_overlaps(rx, rb))
        {
            flint_printf("FAIL: overlap b\n\n");
            flint_printf("a = "); acb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); acb_printn(b, 50, 0); flint_printf("\n\n");
            flint_printf("c = "); acb_printn(c, 50, 0); flint_printf("\n\n");
            flint_printf("x = "); acb_printn(x, 50, ARB_STR_MORE); flint_printf("\n\n");
            flint_printf("rb = "); acb_printn(rb, 50, ARB_STR_MORE); flint_printf("\n\n");
            flint_printf("rx = "); acb_printn(rx, 50, ARB_STR_MORE); flint_printf("\n\n");
            flint_abort();
        }

        if (!acb_overlaps(rx, rc))
        {
            flint_printf("FAIL: overlap c\n\n");
            flint_printf("a = "); acb_printn(a, 50, 0); flint_printf("\n\n");
            flint_printf("b = "); acb_printn(b, 50, 0); flint_printf("\n\n");
            flint_printf("c = "); acb_printn(c, 50, 0); flint_printf("\n\n");
            flint_printf("x = "); acb_printn(x, 50, 0); flint_printf("\n\n");
            flint_printf("rc = "); acb_printn(rc, 50, 0); flint_printf("\n\n");
            flint_printf("rx = "); acb_printn(rx, 50, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
        acb_clear(x);
        acb_clear(ra);
        acb_clear(rb);
        acb_clear(rc);
        acb_clear(rx);
    }

    /* check (a^(-1/2))^(-2) = a */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c;
        slong prec;

        acb_init(a);
        acb_init(b);
        acb_init(c);

        acb_randtest(a, state, 1 + n_randint(state, 2000), 10);
        acb_randtest(b, state, 1 + n_randint(state, 2000), 10);

        prec = 2 + n_randint(state, 2000);

        acb_rsqrt(b, a, prec);
        acb_inv(c, b, prec);
        acb_mul(c, c, c, prec);

        if (!acb_contains(c, a))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        acb_log(c, a, FLINT_MIN(prec, 200));
        acb_mul_2exp_si(c, c, -1);
        acb_neg(c, c);
        acb_exp(c, c, FLINT_MIN(prec, 200));

        if (!acb_overlaps(c, b))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_printf("c = "); acb_print(c); flint_printf("\n\n");
            flint_abort();
        }

        acb_rsqrt(a, a, prec);
        if (!acb_equal(a, b))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("a = "); acb_print(a); flint_printf("\n\n");
            flint_printf("b = "); acb_print(b); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
