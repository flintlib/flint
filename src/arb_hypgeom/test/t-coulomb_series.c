/*
    Copyright (C) 2019 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "arb_poly.h"
#include "arb_hypgeom.h"

TEST_FUNCTION_START(arb_hypgeom_coulomb_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        arb_poly_t F, G, F2, G2, z, w, t, u;
        arb_t c, l, eta, z0;
        slong n1, n2, prec1, prec2;
        unsigned int mask;

        arb_poly_init(F); arb_poly_init(G);
        arb_poly_init(F2); arb_poly_init(G2);
        arb_poly_init(z); arb_poly_init(w);
        arb_poly_init(t);
        arb_poly_init(u);
        arb_init(c); arb_init(l); arb_init(eta); arb_init(z0);

        prec1 = 2 + n_randint(state, 200);
        prec2 = 2 + n_randint(state, 200);

        n1 = n_randint(state, 8);
        n2 = n_randint(state, 8);

        arb_poly_randtest(F, state, 10, prec1, 10);
        arb_poly_randtest(G, state, 10, prec1, 10);

        arb_randtest(l, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        arb_randtest(eta, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        arb_poly_randtest(z, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);

        arb_hypgeom_coulomb_series(F, G, l, eta, z, n1, prec1);

        /* F' G - F G' = 1 */
        arb_poly_derivative(t, F, prec1);
        arb_poly_set(u, G);
        arb_poly_mullow(w, t, u, FLINT_MAX(n1 - 1, 0), prec1);
        arb_poly_derivative(u, u, prec1);
        arb_poly_mullow(t, F, u, FLINT_MAX(n1 - 1, 0), prec1);
        arb_poly_sub(w, w, t, prec1);

        arb_poly_derivative(t, z, prec1);
        arb_poly_truncate(t, FLINT_MAX(n1 - 1, 0));

        /* hack: work around mullow(nan, 0) = 0 */
        arb_poly_get_coeff_arb(z0, z, 0);

        if (!arb_contains_zero(z0) && !arb_poly_overlaps(w, t))
        {
            flint_printf("FAIL: wronskian, n1 = %wd\n\n", n1);
            flint_printf("l = "); arb_printd(l, 30); flint_printf("\n\n");
            flint_printf("eta = "); arb_printd(eta, 30); flint_printf("\n\n");
            flint_printf("z = "); arb_poly_printd(z, 30); flint_printf("\n\n");
            flint_printf("F = "); arb_poly_printd(F, 30); flint_printf("\n\n");
            flint_printf("G = "); arb_poly_printd(G, 30); flint_printf("\n\n");
            flint_printf("w = "); arb_poly_printd(w, 30); flint_printf("\n\n");
            flint_printf("t = "); arb_poly_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        mask = n_randlimb(state);

        arb_poly_set(G2, z); /* for aliasing */

        if (n_randint(state, 2))
        {
            arb_poly_randtest(u, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);
            arb_poly_add(G2, G2, u, prec2);
            arb_poly_sub(G2, G2, u, prec2);
        }

        arb_hypgeom_coulomb_series((mask & 1) ? F2 : NULL,
                                   (mask & 2) ? G2 : NULL, l, eta, G2, n2, prec2);

        arb_poly_truncate(F, FLINT_MIN(n1, n2));
        arb_poly_truncate(G, FLINT_MIN(n1, n2));
        arb_poly_truncate(F2, FLINT_MIN(n1, n2));
        arb_poly_truncate(G2, FLINT_MIN(n1, n2));

        if (((mask & 1) && (!arb_poly_overlaps(F, F2))) ||
            ((mask & 2) && (!arb_poly_overlaps(G, G2))))
        {
            flint_printf("FAIL: consistency (mask)\n\n");
            flint_printf("mask = %u\n\n", mask);
            flint_printf("len1 = %wd, len2 = %wd\n\n", n1, n2);
            flint_printf("l = "); arb_printd(l, 30); flint_printf("\n\n");
            flint_printf("eta = "); arb_printd(eta, 30); flint_printf("\n\n");
            flint_printf("z = "); arb_poly_printd(z, 30); flint_printf("\n\n");
            flint_printf("F = "); arb_poly_printd(F, 30); flint_printf("\n\n");
            flint_printf("F2 = "); arb_poly_printd(F2, 30); flint_printf("\n\n");
            flint_printf("G = "); arb_poly_printd(G, 30); flint_printf("\n\n");
            flint_printf("G2 = "); arb_poly_printd(G2, 30); flint_printf("\n\n");
            flint_abort();
        }

        arb_poly_clear(F); arb_poly_clear(G);
        arb_poly_clear(F2); arb_poly_clear(G2);
        arb_poly_clear(z); arb_poly_clear(w);
        arb_poly_clear(t);
        arb_poly_clear(u);
        arb_clear(c); arb_clear(l); arb_clear(eta); arb_clear(z0);
    }

    TEST_FUNCTION_END(state);
}
