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
#include "acb_poly.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_coulomb_series, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t F, G, Hpos, Hneg, F2, G2, Hpos2, Hneg2, z, w, t, u;
        acb_t c, l, eta, z0;
        slong n1, n2, prec1, prec2;
        unsigned int mask;

        acb_poly_init(F); acb_poly_init(G);
        acb_poly_init(Hpos); acb_poly_init(Hneg);
        acb_poly_init(F2); acb_poly_init(G2);
        acb_poly_init(Hpos2); acb_poly_init(Hneg2);
        acb_poly_init(z); acb_poly_init(w);
        acb_poly_init(t);
        acb_poly_init(u);
        acb_init(c); acb_init(l); acb_init(eta); acb_init(z0);

        prec1 = 2 + n_randint(state, 200);
        prec2 = 2 + n_randint(state, 200);

        n1 = n_randint(state, 8);
        n2 = n_randint(state, 8);

        acb_poly_randtest(F, state, 10, prec1, 10);
        acb_poly_randtest(G, state, 10, prec1, 10);
        acb_poly_randtest(Hpos, state, 10, prec1, 10);
        acb_poly_randtest(Hneg, state, 10, prec1, 10);

        acb_randtest_param(l, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        acb_randtest_param(eta, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        acb_poly_randtest(z, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);

        acb_hypgeom_coulomb_series(F, G, Hpos, Hneg, l, eta, z, n1, prec1);

        acb_poly_derivative(t, F, prec1);

        if (n_randint(state, 2))
            acb_poly_set(u, G);
        else if (n_randint(state, 2))
            acb_poly_set(u, Hpos);
        else
            acb_poly_set(u, Hneg);

        /* F' G - F G' = 1 */
        acb_poly_mullow(w, t, u, FLINT_MAX(n1 - 1, 0), prec1);
        acb_poly_derivative(u, u, prec1);
        acb_poly_mullow(t, F, u, FLINT_MAX(n1 - 1, 0), prec1);
        acb_poly_sub(w, w, t, prec1);

        acb_poly_derivative(t, z, prec1);
        acb_poly_truncate(t, FLINT_MAX(n1 - 1, 0));

        /* hack: work around mullow(nan, 0) = 0 */
        acb_poly_get_coeff_acb(z0, z, 0);

        if (!acb_contains_zero(z0) && !acb_poly_overlaps(w, t))
        {
            flint_printf("FAIL: wronskian, n1 = %wd\n\n", n1);
            flint_printf("l = "); acb_printd(l, 30); flint_printf("\n\n");
            flint_printf("eta = "); acb_printd(eta, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_poly_printd(z, 30); flint_printf("\n\n");
            flint_printf("F = "); acb_poly_printd(F, 30); flint_printf("\n\n");
            flint_printf("G = "); acb_poly_printd(G, 30); flint_printf("\n\n");
            flint_printf("Hpos = "); acb_poly_printd(Hpos, 30); flint_printf("\n\n");
            flint_printf("Hneg = "); acb_poly_printd(Hneg, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_poly_printd(w, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_poly_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        mask = n_randlimb(state);

        acb_poly_set(G2, z); /* for aliasing */

        if (n_randint(state, 2))
        {
            acb_poly_randtest(u, state, 1 + n_randint(state, 10), 1 + n_randint(state, 200), 10);
            acb_poly_add(G2, G2, u, prec2);
            acb_poly_sub(G2, G2, u, prec2);
        }

        acb_hypgeom_coulomb_series((mask & 1) ? F2 : NULL,
                                   (mask & 2) ? G2 : NULL,
                                   (mask & 4) ? Hpos2 : NULL,
                                   (mask & 8) ? Hneg2 : NULL, l, eta, G2, n2, prec2);

        acb_poly_truncate(F, FLINT_MIN(n1, n2));
        acb_poly_truncate(G, FLINT_MIN(n1, n2));
        acb_poly_truncate(Hpos, FLINT_MIN(n1, n2));
        acb_poly_truncate(Hneg, FLINT_MIN(n1, n2));
        acb_poly_truncate(F2, FLINT_MIN(n1, n2));
        acb_poly_truncate(G2, FLINT_MIN(n1, n2));
        acb_poly_truncate(Hpos2, FLINT_MIN(n1, n2));
        acb_poly_truncate(Hneg2, FLINT_MIN(n1, n2));

        if (((mask & 1) && (!acb_poly_overlaps(F, F2))) ||
            ((mask & 2) && (!acb_poly_overlaps(G, G2))) ||
            ((mask & 4) && (!acb_poly_overlaps(Hpos, Hpos2))) ||
            ((mask & 8) && (!acb_poly_overlaps(Hneg, Hneg2))))
        {
            flint_printf("FAIL: consistency (mask)\n\n");
            flint_printf("mask = %u\n\n", mask);
            flint_printf("len1 = %wd, len2 = %wd\n\n", n1, n2);
            flint_printf("l = "); acb_printd(l, 30); flint_printf("\n\n");
            flint_printf("eta = "); acb_printd(eta, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_poly_printd(z, 30); flint_printf("\n\n");
            flint_printf("F = "); acb_poly_printd(F, 30); flint_printf("\n\n");
            flint_printf("F2 = "); acb_poly_printd(F2, 30); flint_printf("\n\n");
            flint_printf("G = "); acb_poly_printd(G, 30); flint_printf("\n\n");
            flint_printf("G2 = "); acb_poly_printd(G2, 30); flint_printf("\n\n");
            flint_printf("Hpos = "); acb_poly_printd(Hpos, 30); flint_printf("\n\n");
            flint_printf("Hpos2 = "); acb_poly_printd(Hpos2, 30); flint_printf("\n\n");
            flint_printf("Hneg = "); acb_poly_printd(Hneg, 30); flint_printf("\n\n");
            flint_printf("Hneg2 = "); acb_poly_printd(Hneg2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(F); acb_poly_clear(G);
        acb_poly_clear(Hpos); acb_poly_clear(Hneg);
        acb_poly_clear(F2); acb_poly_clear(G2);
        acb_poly_clear(Hpos2); acb_poly_clear(Hneg2);
        acb_poly_clear(z); acb_poly_clear(w);
        acb_poly_clear(t);
        acb_poly_clear(u);
        acb_clear(c); acb_clear(l); acb_clear(eta); acb_clear(z0);
    }

    TEST_FUNCTION_END(state);
}
