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

TEST_FUNCTION_START(arb_hypgeom_erf, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);

        arb_init(a);
        arb_init(b);
        arb_init(c);

        arb_randtest_special(a, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(c, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        switch (n_randint(state, 2))
        {
            case 0:
                if (!arb_hypgeom_erf_bb(b, a, 0, prec1))
                    arb_hypgeom_erf(b, a, prec1);
                break;
            default:
                arb_hypgeom_erf(b, a, prec1);
                break;
        }

        switch (n_randint(state, 2))
        {
            case 0:
                if (!arb_hypgeom_erf_bb(c, a, 0, prec2))
                    arb_hypgeom_erf(c, a, prec2);
                break;
            default:
                arb_hypgeom_erf(c, a, prec2);
                break;
        }

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 30); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        slong prec1, prec2;

        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);

        arb_init(a);
        arb_init(b);
        arb_init(c);

        arb_randtest_special(a, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(c, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        switch (n_randint(state, 2))
        {
            case 0:
                if (!arb_hypgeom_erf_bb(b, a, 1, prec1))
                    arb_hypgeom_erfc(b, a, prec1);
                break;
            default:
                arb_hypgeom_erfc(b, a, prec1);
                break;
        }

        switch (n_randint(state, 2))
        {
            case 0:
                if (!arb_hypgeom_erf_bb(c, a, 1, prec2))
                    arb_hypgeom_erfc(c, a, prec2);
                break;
            default:
                arb_hypgeom_erfc(c, a, prec2);
                break;
        }

        if (!arb_overlaps(b, c))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("a = "); arb_printd(a, 30); flint_printf("\n\n");
            flint_printf("b = "); arb_printd(b, 30); flint_printf("\n\n");
            flint_printf("c = "); arb_printd(c, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

#if 0
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        slong prec1;

        prec1 = 2 + n_randint(state, 5000);

        arb_init(a);
        arb_init(b);
        arb_init(c);

        arb_randtest_special(a, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(c, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        arb_hypgeom_erfc(b, a, prec1 + 100);
        arb_set_round(b, b, prec1);
        arb_div_ui(b, b, 3, prec1);
        arb_mul_ui(b, b, 3, prec1);

        if (n_randint(state, 2))
        {
            if (!arb_hypgeom_erf_bb(c, a, 1, prec1))
                arb_hypgeom_erfc(c, a, prec1);
        }
        else
            arb_hypgeom_erfc(c, a, prec1);

        if (arb_is_finite(b) && (arb_rel_accuracy_bits(c) < arb_rel_accuracy_bits(b) - 4.0))
        {
            flint_printf("ACCURACY (erfc)\n\n");
            flint_printf("prec = %wd\n\n", prec1);
            flint_printf("a = "); arb_printd(a, 200); flint_printf("\n");
            flint_printf("b = "); arb_printd(b, 200); flint_printf("\n");
            flint_printf("c = "); arb_printd(c, 200); flint_printf("\n\n");
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

#endif

#if 0
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_t a, b, c;
        slong prec1;

        prec1 = 2 + n_randint(state, 5000);

        arb_init(a);
        arb_init(b);
        arb_init(c);

        arb_randtest_special(a, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        arb_randtest_special(c, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        arb_hypgeom_erf(b, a, prec1 + 100);
        arb_set_round(b, b, prec1);

        arb_div_ui(b, b, 3, prec1);
        arb_mul_ui(b, b, 3, prec1);

        if (n_randint(state, 2))
        {
            if (!arb_hypgeom_erf_bb(c, a, 0, prec1))
                arb_hypgeom_erf(c, a, prec1);
        }
        else
            arb_hypgeom_erf(c, a, prec1);

        if (arb_is_finite(b) && (arb_rel_accuracy_bits(c) < arb_rel_accuracy_bits(b) - 4.0))
        {
            flint_printf("ACCURACY (erf)\n\n");
            flint_printf("prec = %wd\n\n", prec1);
            flint_printf("a = "); arb_printd(a, 200); flint_printf("\n");
            flint_printf("b = "); arb_printd(b, 200); flint_printf("\n");
            flint_printf("c = "); arb_printd(c, 200); flint_printf("\n\n");
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
    }

#endif

    TEST_FUNCTION_END(state);
}
