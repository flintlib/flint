/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_hypgeom.h"

TEST_FUNCTION_START(acb_hypgeom_bessel_k, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t nu0, nu1, nu2, z, w0, w1, w2, t, u;
        slong prec0, prec1, prec2;
        int scaled;

        acb_init(nu0);
        acb_init(nu1);
        acb_init(nu2);
        acb_init(z);
        acb_init(w0);
        acb_init(w1);
        acb_init(w2);
        acb_init(t);
        acb_init(u);

        prec0 = 2 + n_randint(state, 700);
        prec1 = 2 + n_randint(state, 700);
        prec2 = 2 + n_randint(state, 700);

        scaled = n_randint(state, 2);
        acb_randtest_param(nu0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w2, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        if (n_randint(state, 4) == 0)
            arb_zero(acb_imagref(z));

        if (n_randint(state, 2) && prec0 <= 200 && prec1 <= 200 && prec2 <= 200)
        {
            arb_zero(acb_imagref(nu0));
            arb_zero(acb_imagref(z));
            arb_abs(acb_realref(nu0), acb_realref(nu0));
        }

        acb_sub_ui(nu1, nu0, 1, prec0);
        acb_sub_ui(nu2, nu0, 2, prec0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_k_asymp(w0, nu0, z, scaled, prec0);
                break;
            case 1:
                acb_hypgeom_bessel_k_0f1(w0, nu0, z, scaled, prec0);
                break;
            default:
                if (scaled)
                    acb_hypgeom_bessel_k_scaled(w0, nu0, z, prec0);
                else
                    acb_hypgeom_bessel_k(w0, nu0, z, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_k_asymp(w1, nu0, z, scaled, prec1);
                break;
            case 1:
                acb_hypgeom_bessel_k_0f1(w1, nu0, z, scaled, prec1);
                break;
            default:
                if (scaled)
                    acb_hypgeom_bessel_k_scaled(w1, nu0, z, prec1);
                else
                    acb_hypgeom_bessel_k(w1, nu0, z, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("scaled = %d\n\n", scaled);
            flint_printf("nu = "); acb_printd(nu0, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_abort();
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_k_asymp(w1, nu1, z, scaled, prec1);
                break;
            case 1:
                acb_hypgeom_bessel_k_0f1(w1, nu1, z, scaled, prec1);
                break;
            default:
                if (scaled)
                    acb_hypgeom_bessel_k_scaled(w1, nu1, z, prec1);
                else
                    acb_hypgeom_bessel_k(w1, nu1, z, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_k_asymp(w2, nu2, z, scaled, prec2);
                break;
            case 1:
                acb_hypgeom_bessel_k_0f1(w2, nu2, z, scaled, prec2);
                break;
            default:
                if (scaled)
                {
                    acb_hypgeom_bessel_k(w2, nu2, z, prec2);
                    acb_exp(t, z, prec2);
                    acb_mul(w2, w2, t, prec2);
                }
                else
                {
                    acb_hypgeom_bessel_k_scaled(w2, nu2, z, prec2);
                    acb_neg(t, z);
                    acb_exp(t, t, prec2);
                    acb_mul(w2, w2, t, prec2);
                }
        }

        acb_mul(t, w1, nu1, prec0);
        acb_mul_2exp_si(t, t, 1);
        acb_addmul(t, w2, z, prec0);
        acb_submul(t, w0, z, prec0);

        if (!acb_contains_zero(t))
        {
            flint_printf("FAIL: contiguous relation\n\n");
            flint_printf("nu = "); acb_printd(nu0, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(nu0);
        acb_clear(nu1);
        acb_clear(nu2);
        acb_clear(z);
        acb_clear(w0);
        acb_clear(w1);
        acb_clear(w2);
        acb_clear(t);
        acb_clear(u);
    }

    TEST_FUNCTION_END(state);
}
