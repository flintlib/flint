/*
    Copyright (C) 2015 Fredrik Johansson

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

TEST_FUNCTION_START(acb_hypgeom_airy_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t ai, aip, bi, bip, ai2, aip2, bi2, bip2, z, w, t;
        acb_t c;
        slong n1, n2, prec1, prec2;
        unsigned int mask;

        acb_poly_init(ai); acb_poly_init(aip);
        acb_poly_init(bi); acb_poly_init(bip);
        acb_poly_init(ai2); acb_poly_init(aip2);
        acb_poly_init(bi2); acb_poly_init(bip2);
        acb_poly_init(z); acb_poly_init(w);
        acb_poly_init(t); acb_init(c);

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);

        n1 = n_randint(state, 6);
        n2 = n_randint(state, 6);

        acb_poly_randtest(ai, state, 10, prec1, 10);
        acb_poly_randtest(aip, state, 10, prec1, 10);
        acb_poly_randtest(bi, state, 10, prec1, 10);
        acb_poly_randtest(bip, state, 10, prec1, 10);
        acb_poly_randtest(z, state, 1 + n_randint(state, 10), prec1, 10);

        acb_hypgeom_airy_series(ai, aip, bi, bip, z, n1, prec1);

        acb_poly_mullow(w, ai, bip, n1, prec1);
        acb_poly_mullow(t, bi, aip, n1, prec1);
        acb_poly_sub(w, w, t, prec1);

        acb_const_pi(c, prec1);
        acb_inv(c, c, prec1);
        acb_poly_set_acb(t, c);
        acb_poly_truncate(t, n1);

        if (!acb_poly_overlaps(w, t))
        {
            flint_printf("FAIL: wronskian\n\n");
            flint_printf("z = "); acb_poly_printd(z, 30); flint_printf("\n\n");
            flint_printf("ai = "); acb_poly_printd(ai, 30); flint_printf("\n\n");
            flint_printf("aip = "); acb_poly_printd(aip, 30); flint_printf("\n\n");
            flint_printf("bi = "); acb_poly_printd(bi, 30); flint_printf("\n\n");
            flint_printf("bip = "); acb_poly_printd(bip, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_poly_printd(w, 30); flint_printf("\n\n");
            flint_abort();
        }

        mask = n_randlimb(state);

        acb_hypgeom_airy_series((mask & 1) ? ai2 : NULL,
                                (mask & 2) ? aip2 : NULL,
                                (mask & 4) ? bi2 : NULL,
                                (mask & 8) ? bip2 : NULL, z, n2, prec2);

        acb_poly_truncate(ai, FLINT_MIN(n1, n2));
        acb_poly_truncate(aip, FLINT_MIN(n1, n2));
        acb_poly_truncate(bi, FLINT_MIN(n1, n2));
        acb_poly_truncate(bip, FLINT_MIN(n1, n2));
        acb_poly_truncate(ai2, FLINT_MIN(n1, n2));
        acb_poly_truncate(aip2, FLINT_MIN(n1, n2));
        acb_poly_truncate(bi2, FLINT_MIN(n1, n2));
        acb_poly_truncate(bip2, FLINT_MIN(n1, n2));

        if (((mask & 1) && (!acb_poly_overlaps(ai, ai2))) ||
            ((mask & 2) && (!acb_poly_overlaps(aip, aip2))) ||
            ((mask & 4) && (!acb_poly_overlaps(bi, bi2))) ||
            ((mask & 8) && (!acb_poly_overlaps(bip, bip2))))
        {
            flint_printf("FAIL: consistency (mask)\n\n");
            flint_printf("mask = %u\n\n", mask);
            flint_printf("len1 = %wd, len2 = %wd\n\n", n1, n2);
            flint_printf("z = "); acb_poly_printd(z, 30); flint_printf("\n\n");
            flint_printf("ai = "); acb_poly_printd(ai, 30); flint_printf("\n\n");
            flint_printf("ai2 = "); acb_poly_printd(ai2, 30); flint_printf("\n\n");
            flint_printf("aip = "); acb_poly_printd(aip, 30); flint_printf("\n\n");
            flint_printf("aip2 = "); acb_poly_printd(aip2, 30); flint_printf("\n\n");
            flint_printf("bi = "); acb_poly_printd(bi, 30); flint_printf("\n\n");
            flint_printf("bi2 = "); acb_poly_printd(bi2, 30); flint_printf("\n\n");
            flint_printf("bip = "); acb_poly_printd(bip, 30); flint_printf("\n\n");
            flint_printf("bip2 = "); acb_poly_printd(bip2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(ai); acb_poly_clear(aip);
        acb_poly_clear(bi); acb_poly_clear(bip);
        acb_poly_clear(ai2); acb_poly_clear(aip2);
        acb_poly_clear(bi2); acb_poly_clear(bip2);
        acb_poly_clear(z); acb_poly_clear(w);
        acb_poly_clear(t); acb_clear(c);
    }

    TEST_FUNCTION_END(state);
}
