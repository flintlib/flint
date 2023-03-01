/*
    Copyright (C) 2020 D.H.J. Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("bessel_y....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        arb_t x, y, nu;
        slong prec;

        arb_init(x);
        arb_init(y);
        arb_init(nu);

        prec = 2 + n_randint(state, 200);

        arb_randtest(nu, state, 1 + n_randint(state, 200), 1 + n_randint(state, 100));

        arb_pos_inf(x);
        arb_hypgeom_bessel_y(y, nu, x, prec);
        if (arb_is_finite(nu) && !arb_is_zero(y))
        {
            flint_printf("FAIL: positive infinity\n\n");
            flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
            flint_printf("nu = "); arb_printd(nu, 50); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 50); flint_printf("\n\n");
            flint_abort();
        }

        arb_neg_inf(x);
        arb_hypgeom_bessel_y(y, nu, x, prec);
        if (arb_is_finite(nu) && !arb_is_zero(y))
        {
            flint_printf("FAIL: negative infinity\n\n");
            flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
            flint_printf("nu = "); arb_printd(nu, 50); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 50); flint_printf("\n\n");
            flint_abort();
        }

        arb_indeterminate(x);
        arb_hypgeom_bessel_y(y, nu, x, prec);
        if (!arf_is_nan(arb_midref(y)))
        {
            flint_printf("FAIL: indeterminate\n\n");
            flint_printf("x = "); arb_printd(x, 50); flint_printf("\n\n");
            flint_printf("nu = "); arb_printd(nu, 50); flint_printf("\n\n");
            flint_printf("y = "); arb_printd(y, 50); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(nu);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
