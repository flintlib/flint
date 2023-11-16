/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_poly.h"
#include "acb_calc.h"

/* sin(x), defined in t-cauchy_bound.c and t-integrate_taylor.c */
#ifndef sin_x
#define sin_x sin_x
int
sin_x(acb_ptr out, const acb_t inp, void * params, slong order, slong prec)
{
    int xlen = FLINT_MIN(2, order);

    acb_set(out, inp);
    if (xlen > 1)
        acb_one(out + 1);

    _acb_poly_sin_series(out, out, xlen, order, prec);
    return 0;
}
#endif

TEST_FUNCTION_START(acb_calc_integrate_taylor, state)
{
    slong iter;

    for (iter = 0; iter < 150 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t ans, res, a, b;
        arf_t inr, outr;
        double t;
        slong goal, prec;

        acb_init(ans);
        acb_init(res);
        acb_init(a);
        acb_init(b);
        arf_init(inr);
        arf_init(outr);

        goal = 2 + n_randint(state, 300);
        prec = 2 + n_randint(state, 300);

        acb_randtest(a, state, 1 + n_randint(state, 200), 2);
        acb_randtest(b, state, 1 + n_randint(state, 200), 2);

        acb_cos(ans, a, prec);
        acb_cos(res, b, prec);
        acb_sub(ans, ans, res, prec);

        t = (1 + n_randint(state, 20)) / 10.0;
        arf_set_d(inr, t);
        arf_set_d(outr, t + (1 + n_randint(state, 20)) / 5.0);

        acb_calc_integrate_taylor(res, sin_x, NULL,
            a, b, inr, outr, goal, prec);

        if (!acb_overlaps(res, ans))
        {
            flint_printf("FAIL! (iter = %wd)\n", iter);
            flint_printf("prec = %wd, goal = %wd\n", prec, goal);
            flint_printf("inr = "); arf_printd(inr, 15); flint_printf("\n");
            flint_printf("outr = "); arf_printd(outr, 15); flint_printf("\n");
            flint_printf("a = "); acb_printd(a, 15); flint_printf("\n");
            flint_printf("b = "); acb_printd(b, 15); flint_printf("\n");
            flint_printf("res = "); acb_printd(res, 15); flint_printf("\n\n");
            flint_printf("ans = "); acb_printd(ans, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(ans);
        acb_clear(res);
        acb_clear(a);
        acb_clear(b);
        arf_clear(inr);
        arf_clear(outr);
    }

    TEST_FUNCTION_END(state);
}
