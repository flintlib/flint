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

void
acb_hypgeom_u_asymp_proper(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_pow(t, z, a, prec);
    acb_hypgeom_u_asymp(res, a, b, z, -1, prec);
    acb_div(res, res, t, prec);
    acb_clear(t);
}

TEST_FUNCTION_START(acb_hypgeom_u, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t a0, a1, a2, b, z, w0, w1, w2, t, u;
        slong prec0, prec1, prec2;

        acb_init(a0);
        acb_init(a1);
        acb_init(a2);
        acb_init(b);
        acb_init(z);
        acb_init(w0);
        acb_init(w1);
        acb_init(w2);
        acb_init(t);
        acb_init(u);

        prec0 = 2 + n_randint(state, 700);
        prec1 = 2 + n_randint(state, 700);
        prec2 = 2 + n_randint(state, 700);

        acb_randtest_param(a0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest_param(b, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w2, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        if (prec0 <= 300 && prec1 <= 300 && prec2 <= 300 && n_randint(state, 2))
        {
            arb_zero(acb_imagref(a0));
            arb_zero(acb_imagref(b));
            arb_zero(acb_imagref(z));
        }

        acb_add_ui(a1, a0, 1, prec0);
        acb_add_ui(a2, a0, 2, prec0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_u_asymp_proper(w0, a0, b, z, prec0);
                break;
            case 1:
                acb_hypgeom_u_1f1(w0, a0, b, z, prec0);
                break;
            default:
                acb_hypgeom_u(w0, a0, b, z, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_u_asymp_proper(w1, a0, b, z, prec1);
                break;
            case 1:
                acb_hypgeom_u_1f1(w1, a0, b, z, prec1);
                break;
            default:
                acb_hypgeom_u(w1, a0, b, z, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("a = "); acb_printd(a0, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_abort();
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_u_asymp_proper(w1, a1, b, z, prec1);
                break;
            case 1:
                acb_hypgeom_u_1f1(w1, a1, b, z, prec1);
                break;
            default:
                acb_hypgeom_u(w1, a1, b, z, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_u_asymp_proper(w2, a2, b, z, prec2);
                break;
            case 1:
                acb_hypgeom_u_1f1(w2, a2, b, z, prec2);
                break;
            default:
                acb_hypgeom_u(w2, a2, b, z, prec2);
        }

        acb_set(t, w0);

        acb_mul_2exp_si(u, a0, 1);
        acb_sub(u, u, b, prec0);
        acb_add(u, u, z, prec0);
        acb_add_ui(u, u, 2, prec0);

        acb_submul(t, w1, u, prec0);

        acb_sub(u, a2, b, prec0);
        acb_mul(u, u, a1, prec0);
        acb_addmul(t, w2, u, prec0);

        if (!acb_contains_zero(t))
        {
            flint_printf("FAIL: contiguous relation\n\n");
            flint_printf("a = "); acb_printd(a0, 30); flint_printf("\n\n");
            flint_printf("b = "); acb_printd(b, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("w2 = "); acb_printd(w2, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(a0);
        acb_clear(a1);
        acb_clear(a2);
        acb_clear(b);
        acb_clear(z);
        acb_clear(w0);
        acb_clear(w1);
        acb_clear(w2);
        acb_clear(t);
        acb_clear(u);
    }

    TEST_FUNCTION_END(state);
}
