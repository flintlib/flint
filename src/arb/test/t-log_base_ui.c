/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "arb.h"

TEST_FUNCTION_START(arb_log_base_ui, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, bx, lbx, lbx1;
        ulong n, b, b2, r;
        slong prec;

        arb_init(x);
        arb_init(bx);
        arb_init(lbx);
        arb_init(lbx1);

        prec = 2 + n_randint(state, 500);

        b = n_randtest(state);
        b = n_pow(b, 1 + n_randint(state, 4));

        if (n_randint(state, 2))
            arb_set_ui(x, n_randtest(state));
        else
            arb_randtest(x, state, 2 + n_randint(state, 500), 1 + n_randint(state, 200));

        arb_randtest(bx, state, 2 + n_randint(state, 500), 1 + n_randint(state, 200));
        arb_randtest(lbx, state, 2 + n_randint(state, 500), 1 + n_randint(state, 200));

        arb_set_ui(bx, b);
        arb_pow(bx, bx, x, 2 + n_randint(state, 500));
        arb_log_base_ui(lbx, bx, b, prec);

        /* test log_b(b^x) = x */
        if (!arb_overlaps(lbx, x))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("b = %wu\n\n", b);
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
            flint_printf("bx = "); arb_printd(bx, 30); flint_printf("\n\n");
            flint_printf("lbx = "); arb_printd(lbx, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_log_base_ui(bx, bx, b, prec);

        if (!arb_overlaps(bx, lbx))
        {
            flint_printf("FAIL: aliasing\n\n");
            flint_printf("b = %wu\n\n", b);
            flint_printf("x = "); arb_printd(x, 30); flint_printf("\n\n");
            flint_printf("bx = "); arb_printd(bx, 30); flint_printf("\n\n");
            flint_printf("lbx = "); arb_printd(lbx, 30); flint_printf("\n\n");
            flint_abort();
        }

        /* test exact computation of log_{b^(2^r)}(b^n)*/
        n = n_randint(state, 100);
        b = 2 + n_randint(state, 10);
        r = n_randint(state, 4);
        b2 = n_pow(b, 1 << r);

        arb_set_ui(lbx1, n);
        arb_mul_2exp_si(lbx1, lbx1, -r);

        for (prec = 20; ; prec *= 2)
        {
            arb_set_ui(bx, b);
            arb_pow_ui(bx, bx, n, prec);
            arb_log_base_ui(lbx, bx, b2, prec);

            if (!arb_contains(lbx, lbx1) || prec > 10000)
            {
                flint_printf("FAIL: containment or exactness\n\n");
                flint_printf("b = %wu\n\n", b);
                flint_printf("b2 = %wu\n\n", b2);
                flint_printf("r = %wu\n\n", r);
                flint_printf("n = %wu\n\n", n);
                flint_printf("bx = "); arb_printd(bx, 30); flint_printf("\n\n");
                flint_printf("lbx = "); arb_printd(lbx, 30); flint_printf("\n\n");
                flint_printf("lbx1 = "); arb_printd(lbx1, 30); flint_printf("\n\n");
                flint_abort();
            }

            if (arb_is_exact(lbx))
                break;
        }

        arb_clear(x);
        arb_clear(bx);
        arb_clear(lbx);
        arb_clear(lbx1);
    }

    TEST_FUNCTION_END(state);
}
