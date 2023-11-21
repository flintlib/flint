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

TEST_FUNCTION_START(arb_div_newton, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c, d;
        slong prec;

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(d);

        prec = 2 + n_randint(state, 200000);

        arb_randtest(a, state, 1 + n_randint(state, 200000), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200000), 10);
        arb_randtest(c, state, 1 + n_randint(state, 200000), 10);

        if (n_randint(state, 2))
        {
            arb_set(c, a);
            arb_div_newton(c, c, b, prec);
        }
        else if (n_randint(state, 2))
        {
            arb_set(c, b);
            arb_div_newton(c, a, c, prec);
        }
        else
        {
            arb_div_newton(c, a, b, prec);
        }

        arb_div(d, a, b, prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("a = "); arb_printd(a, 100); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 100); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 100); flint_printf("\n\n");
            flint_printf("d = "); arb_printd(d, 100); flint_printf("\n\n");
            flint_abort();
        }

        arb_randtest(c, state, 1 + n_randint(state, 100000), 10);

        mag_zero(arb_radref(b));

        if (n_randint(state, 2))
        {
            arb_set(c, a);
            arb_div_arf_newton(c, c, arb_midref(b), prec);
        }
        else if (n_randint(state, 2))
        {
            arb_set(c, b);
            arb_div_arf_newton(c, a, arb_midref(c), prec);
        }
        else
        {
            arb_div_arf_newton(c, a, arb_midref(b), prec);
        }

        arb_div_arf(d, a, arb_midref(b), prec);

        if (!arb_overlaps(c, d))
        {
            flint_printf("FAIL: containment\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("a = "); arb_printd(a, 100); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 100); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 100); flint_printf("\n\n");
            flint_printf("d = "); arb_printd(d, 100); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(d);
    }

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpz_t a, b, c, d, e;

        fmpz_init(a);
        fmpz_init(b);
        fmpz_init(c);
        fmpz_init(d);
        fmpz_init(e);

        fmpz_randtest(a, state, 400000);
        fmpz_randtest_not_zero(b, state, 200000);

        arb_fmpz_divapprox(c, a, b);

        fmpz_fdiv_q(d, a, b);
        fmpz_cdiv_q(e, a, b);

        if (!(fmpz_equal(c, d) || fmpz_equal(c, e)))
        {
            flint_printf("FAIL: fmpz_divapprox\n\n");
            flint_printf("a = "); fmpz_print(a); flint_printf("\n\n");
            flint_printf("b = "); fmpz_print(b); flint_printf("\n\n");
            flint_printf("c = "); fmpz_print(c); flint_printf("\n\n");
            flint_abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_clear(e);
    }

    TEST_FUNCTION_END(state);
}
