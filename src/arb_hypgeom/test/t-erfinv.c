/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_erfinv, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z, w;
        slong prec, prec2;

        arb_init(x);
        arb_init(y);
        arb_init(z);
        arb_init(w);

        prec = 2 + n_randint(state, 500);
        prec2 = 2 + n_randint(state, 500);

        arb_init(x);
        arb_init(y);

        arb_randtest(x, state, prec, 10 + n_randint(state, 300));
        arb_randtest(y, state, prec, 100);

        arb_hypgeom_erfcinv(y, x, prec);
        arb_hypgeom_erfc(z, y, prec + 10);

        if (!arb_overlaps(x, z))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("x = "); arb_printn(x, 100, 0); flint_printf("\n\n");
            flint_printf("y = "); arb_printn(y, 100, 0); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_is_exact(x) && arf_sgn(arb_midref(x)) > 0 &&
            arf_cmp_2exp_si(arb_midref(x), 1) < 0 &&
            arb_rel_accuracy_bits(y) < prec - 2)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("iter = %wd, prec = %wd\n", iter, prec);
            flint_printf("x = "); arb_printn(x, 100, 0); flint_printf("\n\n");
            flint_printf("y = "); arb_printn(y, 100, 0); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_set(z, x);
        mag_zero(arb_radref(z));

        if (n_randint(state, 2))
        {
            arf_set_mag(arb_midref(w), arb_radref(x));

            if (n_randint(state, 2))
                arb_neg(w, w);

            arb_add(z, z, w, prec);
        }

        arb_hypgeom_erfcinv(z, z, prec2);

        if (!arb_overlaps(y, z))
        {
            flint_printf("FAIL: overlap 2\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("x = "); arb_printn(x, 100, 0); flint_printf("\n\n");
            flint_printf("y = "); arb_printn(y, 100, 0); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
        arb_clear(w);
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t x, y, z, w;
        slong prec, prec2;

        arb_init(x);
        arb_init(y);
        arb_init(z);
        arb_init(w);

        prec = 2 + n_randint(state, 500);
        prec2 = 2 + n_randint(state, 500);

        arb_init(x);
        arb_init(y);

        arb_randtest(x, state, prec, 100);
        arb_randtest(y, state, prec, 100);

        arb_hypgeom_erfinv(y, x, prec);
        arb_hypgeom_erf(z, y, prec + 10);

        if (!arb_overlaps(x, z))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("x = "); arb_printn(x, 100, 0); flint_printf("\n\n");
            flint_printf("y = "); arb_printn(y, 100, 0); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        if (arb_is_exact(x) && arf_cmpabs_2exp_si(arb_midref(x), 0) < 0 &&
            arb_rel_accuracy_bits(y) < prec - 2)
        {
            flint_printf("FAIL: accuracy\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("x = "); arb_printn(x, 100, 0); flint_printf("\n\n");
            flint_printf("y = "); arb_printn(y, 100, 0); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_set(z, x);
        mag_zero(arb_radref(z));

        if (n_randint(state, 2))
        {
            arf_set_mag(arb_midref(w), arb_radref(x));

            if (n_randint(state, 2))
                arb_neg(w, w);

            arb_add(z, z, w, prec);
        }

        arb_hypgeom_erfinv(z, z, prec2);

        if (!arb_overlaps(y, z))
        {
            flint_printf("FAIL: overlap 2\n\n");
            flint_printf("prec = %wd\n", prec);
            flint_printf("x = "); arb_printn(x, 100, 0); flint_printf("\n\n");
            flint_printf("y = "); arb_printn(y, 100, 0); flint_printf("\n\n");
            flint_printf("z = "); arb_printn(z, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
        arb_clear(w);
    }

    TEST_FUNCTION_END(state);
}
