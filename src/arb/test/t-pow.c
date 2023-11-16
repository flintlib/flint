/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "double_extras.h"

TEST_FUNCTION_START(arb_pow, state)
{
    slong iter;

    /* check large arguments */
    for (iter = 0; iter < 20000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d, e, f;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);
        arb_init(e);
        arb_init(f);

        arb_randtest(a, state, 1 + n_randint(state, 1000), 200);
        arb_randtest(b, state, 1 + n_randint(state, 1000), 200);

        arb_pow(c, a, b, prec1);
        arb_pow(d, a, b, prec2);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_randtest(c, state, 1 + n_randint(state, 1000), 200);

        /* check a^(b+c) = a^b*a^c */
        arb_add(e, b, c, prec1);
        arb_pow(d, a, e, prec1);

        arb_pow(e, a, b, prec1);
        arb_pow(f, a, c, prec1);
        arb_mul(e, e, f, prec1);

        if (!arb_overlaps(d, e))
        {
            flint_printf("FAIL: functional equation\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_printf("e = "); arb_print(e); flint_printf("\n\n");
            flint_abort();
        }

        arb_pow(c, a, b, prec1);
        arb_set(d, a);
        arb_pow(d, d, b, prec2);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: aliasing 1\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_set(d, b);
        arb_pow(d, a, d, prec2);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: aliasing 2\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
        arb_clear(e);
        arb_clear(f);
    }

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        ulong n, prec;
        double epsilon, delta;
        arb_t a, b, c, d;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        prec = 2 + n_randint(state, 100);

        mag_randtest(arb_radref(a), state, 1 + n_randint(state, 1000));
        arf_set_mag(arb_midref(a), arb_radref(a));

        n = 1 + n_randint(state, 100);
        epsilon = d_randtest(state) / 16; /* 1/32 <= epsilon < 1/16 */
        if (n_randint(state, 2))
            epsilon = -epsilon;
        delta = d_randtest(state) / 4; /* 1/8 <= delta < 1/4 */
        arf_set_d(arb_midref(b), n + epsilon);
        mag_set_d(arb_radref(b), delta);

        arb_pow(c, a, b, prec);
        arb_pow_ui(d, a, n, prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: interval from zero\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("n = %d\n\n", n);
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_set(c, a);
        arb_pow(c, c, b, prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: interval from zero; alias 1\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("n = %d\n\n", n);
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_set(c, b);
        arb_pow(c, a, c, prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: interval from zero; alias 2\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("n = %d\n\n", n);
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("d = "); arb_print(d); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(d);
        arb_clear(c);
        arb_clear(b);
        arb_clear(a);
    }

    TEST_FUNCTION_END(state);
}
