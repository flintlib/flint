/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"

TEST_FUNCTION_START(arb_sqrt_newton, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, c, d;
        slong prec;

        arb_init(a);
        arb_init(c);
        arb_init(d);

        prec = 2 + n_randint(state, 200000);

        arb_randtest(a, state, 1 + n_randint(state, 200000), 10);
        arb_randtest(c, state, 1 + n_randint(state, 200000), 10);

        if (n_randint(state, 2))
        {
            arb_set(c, a);
            arb_sqrt_newton(c, c, prec);
        }
        else
        {
            arb_sqrt_newton(c, a, prec);
        }

        arb_sqrt(d, a, prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("a = "); arb_printd(a, 100); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 100); flint_printf("\n\n");
            flint_printf("d = "); arb_printd(d, 100); flint_printf("\n\n");
            flint_abort();
        }

        arb_randtest(c, state, 1 + n_randint(state, 100000), 10);

        mag_zero(arb_radref(a));

        if (n_randint(state, 2))
        {
            arb_set(c, a);
            arb_sqrt_arf_newton(c, arb_midref(a), prec);
        }
        else
        {
            arb_sqrt_arf_newton(c, arb_midref(a), prec);
        }

        arb_sqrt_arf(d, arb_midref(a), prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("a = "); arb_printd(a, 100); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 100); flint_printf("\n\n");
            flint_printf("d = "); arb_printd(d, 100); flint_printf("\n\n");
            flint_abort();
        }

        if (n_randint(state, 2))
        {
            arb_set(c, a);
            arb_rsqrt_arf_newton(c, arb_midref(c), prec);
        }
        else
        {
            arb_rsqrt_arf_newton(c, arb_midref(a), prec);
        }

/*
        BUG in MPFR 4.2.0

        arb_rsqrt(d, a, prec);
*/
        arb_sqrt(d, a, prec + 10);
        arb_inv(d, d, prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: containment (rsqrt)\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("a = "); arb_printd(a, 100); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 100); flint_printf("\n\n");
            flint_printf("d = "); arb_printd(d, 100); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(c);
        arb_clear(d);
    }

    TEST_FUNCTION_END(state);
}
