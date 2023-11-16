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
#include "acb_modular.h"

TEST_FUNCTION_START(acb_modular_theta_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_poly_t t1, t2, t3, t4, t1b, t2b, t3b, t4b, z;
        acb_t tau;
        slong len1, len2, prec1, prec2;
        unsigned int mask;

        acb_poly_init(t1); acb_poly_init(t1b);
        acb_poly_init(t2); acb_poly_init(t2b);
        acb_poly_init(t3); acb_poly_init(t3b);
        acb_poly_init(t4); acb_poly_init(t4b);
        acb_poly_init(z);
        acb_init(tau);

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);

        len1 = n_randint(state, 6);
        len2 = n_randint(state, 6);

        acb_poly_randtest(t1, state, 10, prec1, 10);
        acb_poly_randtest(t2, state, 10, prec1, 10);
        acb_poly_randtest(t3, state, 10, prec1, 10);
        acb_poly_randtest(t4, state, 10, prec1, 10);
        acb_poly_randtest(z, state, 1 + n_randint(state, 10), prec1, 10);
        acb_randtest(tau, state, prec1, 10);

        acb_modular_theta_series(t1, t2, t3, t4, z, tau, len1, prec1);

        mask = n_randlimb(state);

        acb_modular_theta_series((mask & 1) ? t1b : NULL,
                                (mask & 2) ? t2b : NULL,
                                (mask & 4) ? t3b : NULL,
                                (mask & 8) ? t4b : NULL, z, tau, len2, prec2);

        acb_poly_truncate(t1, FLINT_MIN(len1, len2));
        acb_poly_truncate(t1b, FLINT_MIN(len1, len2));
        acb_poly_truncate(t2, FLINT_MIN(len1, len2));
        acb_poly_truncate(t2b, FLINT_MIN(len1, len2));
        acb_poly_truncate(t3, FLINT_MIN(len1, len2));
        acb_poly_truncate(t3b, FLINT_MIN(len1, len2));
        acb_poly_truncate(t4, FLINT_MIN(len1, len2));
        acb_poly_truncate(t4b, FLINT_MIN(len1, len2));

        if (((mask & 1) && (!acb_poly_overlaps(t1, t1b))) ||
            ((mask & 2) && (!acb_poly_overlaps(t2, t2b))) ||
            ((mask & 4) && (!acb_poly_overlaps(t3, t3b))) ||
            ((mask & 8) && (!acb_poly_overlaps(t4, t4b))))
        {
            flint_printf("FAIL: consistency (mask)\n\n");
            flint_printf("mask = %u\n\n", mask);
            flint_printf("len1 = %wd, len2 = %wd\n\n", len1, len2);
            flint_printf("z = "); acb_poly_printd(z, 30); flint_printf("\n\n");
            flint_printf("tau = "); acb_printd(tau, 30); flint_printf("\n\n");
            flint_printf("t1 = "); acb_poly_printd(t1, 30); flint_printf("\n\n");
            flint_printf("t1b = "); acb_poly_printd(t1b, 30); flint_printf("\n\n");
            flint_printf("t2 = "); acb_poly_printd(t2, 30); flint_printf("\n\n");
            flint_printf("t2b = "); acb_poly_printd(t2b, 30); flint_printf("\n\n");
            flint_printf("t3 = "); acb_poly_printd(t3, 30); flint_printf("\n\n");
            flint_printf("t3b = "); acb_poly_printd(t3b, 30); flint_printf("\n\n");
            flint_printf("t4 = "); acb_poly_printd(t4, 30); flint_printf("\n\n");
            flint_printf("t4b = "); acb_poly_printd(t4b, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(t1); acb_poly_clear(t1b);
        acb_poly_clear(t2); acb_poly_clear(t2b);
        acb_poly_clear(t3); acb_poly_clear(t3b);
        acb_poly_clear(t4); acb_poly_clear(t4b);
        acb_poly_clear(z); acb_clear(tau);
    }

    TEST_FUNCTION_END(state);
}
