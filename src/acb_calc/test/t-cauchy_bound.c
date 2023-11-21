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

static const double answers[10] = {
  1.04570093561423094, 2.0358667496686487, 4.82706400405656566,
  11.2347033850581819, 27.1331828778516522, 67.210305990439042,
  168.564351176369878, 427.496180295369909, 1093.56959348921777,
  2815.70144392142227
};

TEST_FUNCTION_START(acb_calc_cauchy_bound, state)
{
    slong iter;

    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t b, radius, ans;
        acb_t x;
        slong r, prec, maxdepth;

        arb_init(b);
        arb_init(radius);
        arb_init(ans);
        acb_init(x);

        acb_set_ui(x, 5);

        r = 1 + n_randint(state, 10);
        arb_set_ui(radius, r);

        prec = 2 + n_randint(state, 100);
        maxdepth = n_randint(state, 10);

        acb_calc_cauchy_bound(b, sin_x, NULL, x, radius, maxdepth, prec);

        arf_set_d(arb_midref(ans), answers[r-1]);
        mag_set_d(arb_radref(ans), 1e-8);

        if (!arb_overlaps(b, ans))
        {
            flint_printf("FAIL\n");
            flint_printf("r = %wd, prec = %wd, maxdepth = %wd\n\n", r, prec, maxdepth);
            arb_printd(b, 15); flint_printf("\n\n");
            arb_printd(ans, 15); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(b);
        arb_clear(radius);
        arb_clear(ans);
        acb_clear(x);
    }

    TEST_FUNCTION_END(state);
}
