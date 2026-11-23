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
randc_1f1(acb_t x, flint_rand_t state, double scale, int real)
{
    acb_set_d_d(x, scale * (n_randint(state, 2001) / 1000.0 - 1.0),
        real ? 0.0 : scale * (n_randint(state, 2001) / 1000.0 - 1.0));
}

TEST_FUNCTION_START(acb_hypgeom_1f1_integration, state)
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

        /* good: Re(a) >= 1, Re(b - a) >= 1 */
        randc_1f1(a, state, scale, real);
        randc_1f1(b, state, scale, real);
        randc_1f1(z, state, 2.0 * scale, n_randint(state, 6) == 0);
        if (good)
        {
            arb_abs(acb_realref(a), acb_realref(a));
            arb_abs(acb_realref(b), acb_realref(b));
            arb_add_ui(acb_realref(a), acb_realref(a), 1, 60);
            arb_add_ui(acb_realref(b), acb_realref(b), 1, 60);
            acb_add(b, b, a, 60);
        }
        if (n_randint(state, 5) == 0)
        {
            arb_floor(acb_realref(a), acb_realref(a), 60);
            arb_zero(acb_imagref(a));
            if (good)
                arb_add_ui(acb_realref(a), acb_realref(a), 1, 60);
        }

        if (n_randint(state, 4) == 0)
            mag_set_ui_2exp_si(arb_radref(acb_realref(z)), 1, -20 - (slong) n_randint(state, 60));

        acb_hypgeom_1f1_integration(r1, a, b, z, flags, prec1);

        if (good)
        {
            nice++;
            computed += acb_is_finite(r1);
        }

        if (acb_is_finite(r1))
        {
            if (n_randint(state, 2))
            {
                _acb_hypgeom_m_nointegration(r2, a, b, z, flags, prec2 + 3 * (slong) scale + 100);
            }
            else
            {
                acb_hypgeom_1f1_integration(r2, a, b, z, flags, prec2);
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

    /* a high precision, where the integrand in the outer cells is below
       the range of doubles relative to its maximum */
    {
        acb_t a, b, z, r1, r2;
        acb_init(a); acb_init(b); acb_init(z); acb_init(r1); acb_init(r2);
        acb_set_d_d(a, 1000.25, 200.0);
        acb_set_d_d(b, 3000.75, -100.0);
        acb_set_d_d(z, 2000.0, 1000.0);
        acb_hypgeom_1f1_integration(r1, a, b, z, 0, 1500);
        acb_hypgeom_1f1_integration(r2, a, b, z, 0, 1200);
        if (!acb_overlaps(r1, r2) || acb_rel_accuracy_bits(r1) < 1470 || acb_rel_accuracy_bits(r2) < 1170)
        {
            flint_printf("FAIL: high precision\n\n");
            flint_printf("accuracy: %wd, %wd\n\n", acb_rel_accuracy_bits(r1), acb_rel_accuracy_bits(r2));
            flint_abort();
        }
        acb_clear(a); acb_clear(b); acb_clear(z); acb_clear(r1); acb_clear(r2);
    }

    TEST_FUNCTION_END(state);
}
