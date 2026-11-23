/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.2

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"
#include "acb_hypgeom/impl.h"

/* a random number with real and imaginary parts in [-scale, scale] */
static void
randc_u(acb_t x, flint_rand_t state, double scale, int real)
{
    acb_set_d_d(x, scale * (n_randint(state, 2001) / 1000.0 - 1.0),
        real ? 0.0 : scale * (n_randint(state, 2001) / 1000.0 - 1.0));
}

TEST_FUNCTION_START(acb_hypgeom_u_integration, state)
{
    slong iter, nice = 0, computed = 0;

    for (iter = 0; iter < 300 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a, b, c, z, r1, r2, t;
        slong prec1, prec2;
        double scale;
        int flags, good, real;

        acb_init(a); acb_init(b); acb_init(c); acb_init(z); acb_init(r1); acb_init(r2); acb_init(t);

        prec1 = 20 + n_randint(state, (n_randint(state, 4) == 0) ? 300 : 100);
        prec2 = prec1 + 20 + n_randint(state, 100);
        scale = exp(log(200.0) * n_randint(state, 1001) / 1000.0);
        flags = n_randint(state, 2);
        real = (n_randint(state, 4) == 0);
        good = (n_randint(state, 5) != 0);

        /* good: Re(a) >= 1, z off the negative axis */
        randc_u(a, state, scale, real);
        randc_u(b, state, 2.0 * scale, real);
        randc_u(z, state, 2.0 * scale, n_randint(state, 6) == 0);
        if (good)
        {
            arb_abs(acb_realref(a), acb_realref(a));
            arb_add_ui(acb_realref(a), acb_realref(a), 1, 60);
            if (arb_is_zero(acb_imagref(z)))
                arb_abs(acb_realref(z), acb_realref(z));
            if (acb_is_zero(z))
                acb_one(z);
        }
        if (n_randint(state, 5) == 0)
        {
            arb_floor(acb_realref(b), acb_realref(b), 60);
            arb_zero(acb_imagref(b));
        }

        if (n_randint(state, 4) == 0)
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -20 - (slong) n_randint(state, 60));

        acb_hypgeom_u_integration(r1, a, b, z, prec1);

        if (good)
        {
            nice++;
            computed += acb_is_finite(r1);
        }

        if (acb_is_finite(r1))
        {
            if (n_randint(state, 2))
            {
                acb_hypgeom_u(r2, a, b, z, prec2 + 3 * (slong) scale + 100);
            }
            else
            {
                acb_hypgeom_u_integration(r2, a, b, z, prec2);
            }

            if (!acb_overlaps(r1, r2))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
                flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
                flint_printf("c = "); acb_printd(c, 30); flint_printf("\n\n");
                flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
                flint_printf("flags = %d, prec1 = %wd, prec2 = %wd\n\n", flags, prec1, prec2);
                flint_printf("r1 = "); acb_printd(r1, 30); flint_printf("\n\n");
                flint_printf("r2 = "); acb_printd(r2, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_clear(a); acb_clear(b); acb_clear(c); acb_clear(z); acb_clear(r1); acb_clear(r2); acb_clear(t);
    }

    if (nice >= 20 && computed < 0.8 * nice)
    {
        flint_printf("FAIL: %wd of %wd computed\n\n", computed, nice);
        flint_abort();
    }

    TEST_FUNCTION_END(state);
}
