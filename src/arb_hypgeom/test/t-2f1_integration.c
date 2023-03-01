/*
    Copyright (C) 2021 Fredrik Johansson

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

    flint_printf("2f1_integration....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b, c, z, r1, r2;
        slong prec1, prec2;
        int regularized;

        prec1 = 2 + n_randint(state, 80);
        prec2 = 2 + n_randint(state, 150);

        arb_init(a);
        arb_init(b);
        arb_init(c);
        arb_init(z);
        arb_init(r1);
        arb_init(r2);

        regularized = n_randint(state, 2);
        arb_randtest_precise(a, state, 1 + n_randint(state, 200), 1 + n_randint(state, 8));
        arb_randtest_precise(b, state, 1 + n_randint(state, 200), 1 + n_randint(state, 8));
        arb_randtest_precise(c, state, 1 + n_randint(state, 200), 1 + n_randint(state, 8));
        arb_randtest_precise(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 8));
        arb_randtest(r1, state, 1 + n_randint(state, 200), 1 + n_randint(state, 5));
        arb_randtest(r2, state, 1 + n_randint(state, 200), 1 + n_randint(state, 5));

        arb_add_ui(a, a, n_randint(state, 100), prec1 + 100);
        arb_add_ui(b, b, n_randint(state, 200), prec1 + 100);
        arb_add_ui(c, c, n_randint(state, 300), prec1 + 100);
        arb_add_ui(z, z, n_randint(state, 100), prec1 + 100);

        arb_hypgeom_2f1_integration(r1, a, b, c, z, regularized, prec1);

        if (arb_is_finite(r1))
        {
            if (n_randint(state, 2))
                arb_hypgeom_2f1(r2, a, b, c, z, regularized, prec2);
            else
                arb_hypgeom_2f1_integration(r2, a, b, c, z, regularized, prec2);

            if (!arb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("a = "); arb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); arb_printd(b, 30); flint_printf("\n\n");
                flint_printf("c = "); arb_printd(c, 30); flint_printf("\n\n");
                flint_printf("z = "); arb_printd(z, 30); flint_printf("\n\n");
                flint_printf("r1 = "); arb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); arb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        if (iter == 0)
        {
            prec1 = 333;

            arb_set_str(b, "2651.1", prec1);
            arb_set_str(a, "50131.3", prec1);
            arb_set_str(c, "72923.5", prec1);
            arb_set_str(z, "-123.5", prec1);

            arb_hypgeom_2f1_integration(r1, a, b, c, z, 0, prec1);
            arb_set_str(r2, "2.334513854528911419992516597096504294432862073292366983166036504547076680083826439493425898608e-5118 +/- 2.73e-5212", prec1);

            if (!arb_overlaps(r1, r2) || arb_rel_accuracy_bits(r1) < arb_rel_accuracy_bits(r2) - 10)
            {
                flint_printf("FAIL: overlap (1)\n\n");
                flint_printf("r1 = "); arb_printd(r1, 100); flint_printf("\n\n");
                flint_printf("r2 = "); arb_printd(r2, 100); flint_printf("\n\n");
                flint_abort();
            }

            arb_set_str(b, "1.1e10", prec1);
            arb_set_str(a, "1.2e11", prec1);
            arb_set_str(c, "1.3e12", prec1);
            arb_set_str(z, "0.7", prec1);

            arb_hypgeom_2f1_integration(r1, a, b, c, z, 0, prec1);
            arb_set_str(r2, "5.86281189034275194553834976777402286036616612793483753071706787203236651175920694151330e+320059600 +/- 1.96e+320059513", prec1);

            if (!arb_overlaps(r1, r2) || arb_rel_accuracy_bits(r1) < arb_rel_accuracy_bits(r2) - 10)
            {
                flint_printf("FAIL: overlap (2)\n\n");
                flint_printf("r1 = "); arb_printd(r1, 100); flint_printf("\n\n");
                flint_printf("r2 = "); arb_printd(r2, 100); flint_printf("\n\n");
                flint_abort();
            }
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);
        arb_clear(z);
        arb_clear(r1);
        arb_clear(r2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
