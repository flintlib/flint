/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "arb.h"

TEST_FUNCTION_START(arb_div, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);
        arb_init(c);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        do {
            arb_randtest(a, state, 1 + n_randint(state, 200), 10);
            arb_randtest(b, state, 1 + n_randint(state, 200), 10);
            arb_randtest(c, state, 1 + n_randint(state, 200), 10);

            arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
            arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        } while (fmpq_is_zero(y));

        arb_div(c, a, b, 2 + n_randint(state, 200));
        fmpq_div(z, x, y);

        if (!arb_contains_fmpq(c, z))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("c = "); arb_print(c); flint_printf("\n\n");
            flint_printf("z = "); fmpq_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and a */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        do {
            arb_randtest(a, state, 1 + n_randint(state, 200), 10);
            arb_randtest(b, state, 1 + n_randint(state, 200), 10);

            arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
            arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        } while (fmpq_is_zero(y));

        arb_div(a, a, b, 2 + n_randint(state, 200));
        fmpq_div(z, x, y);

        if (!arb_contains_fmpq(a, z))
        {
            flint_printf("FAIL: aliasing (c, a)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpq_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and b */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        do {
            arb_randtest(a, state, 1 + n_randint(state, 200), 10);
            arb_randtest(b, state, 1 + n_randint(state, 200), 10);

            arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
            arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        } while (fmpq_is_zero(y));

        arb_div(b, a, b, 2 + n_randint(state, 200));
        fmpq_div(z, x, y);

        if (!arb_contains_fmpq(b, z))
        {
            flint_printf("FAIL: aliasing (c, b)\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("y = "); fmpq_print(y); flint_printf("\n\n");
            flint_printf("z = "); fmpq_print(z); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* test special values */
    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        arb_randtest_special(a, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        arb_randtest_special(b, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));
        arb_randtest_special(c, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));

        arb_div(c, a, b, 2 + n_randint(state, 200));
        arb_mul(d, c, b, 2 + n_randint(state, 200));

        if (!arb_contains(d, a))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 15); flint_printf("\n\n");
            flint_printf("d = "); arb_printd(d, 15); flint_printf("\n\n");
            flint_abort();
        }

        if (arf_is_nan(arb_midref(a)) || arf_is_nan(arb_midref(b)) ||
            arb_contains_zero(b) || (!arb_is_finite(a) && !arb_is_finite(b)))
        {
            if (!arf_is_nan(arb_midref(c)) || !mag_is_inf(arb_radref(c)))
            {
                flint_printf("FAIL: special value 1\n\n");
                flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
                flint_printf("b = "); arb_printd(b, 15); flint_printf("\n\n");
                flint_printf("c = "); arb_printd(c, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (!arb_is_finite(a) && !arb_contains_zero(a) && !arb_contains_zero(b)
            && arb_is_finite(b))
        {
            if (!arf_is_inf(arb_midref(c)) || !mag_is_zero(arb_radref(c)) ||
                arf_sgn(arb_midref(a)) * arf_sgn(arb_midref(b)) != arf_sgn(arb_midref(c)))
            {
                flint_printf("FAIL: special value 2\n\n");
                flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
                flint_printf("b = "); arb_printd(b, 15); flint_printf("\n\n");
                flint_printf("c = "); arb_printd(c, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (!arb_is_finite(a) && !arf_is_nan(arb_midref(a)) &&
            arb_contains_zero(a) &&
            !arb_contains_zero(b) && arb_is_finite(b))
        {
            if (!arf_is_zero(arb_midref(c)) || !mag_is_inf(arb_radref(c)))
            {
                flint_printf("FAIL: special value 3\n\n");
                flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
                flint_printf("b = "); arb_printd(b, 15); flint_printf("\n\n");
                flint_printf("c = "); arb_printd(c, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (arb_is_finite(a) && !arb_is_finite(b) && !arb_contains_zero(b))
        {
            if (!arb_is_zero(c))
            {
                flint_printf("FAIL: special value 4\n\n");
                flint_printf("a = "); arb_printd(a, 15); flint_printf("\n\n");
                flint_printf("b = "); arb_printd(b, 15); flint_printf("\n\n");
                flint_printf("c = "); arb_printd(c, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    TEST_FUNCTION_END(state);
}
