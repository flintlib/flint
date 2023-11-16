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
#include "acb_modular.h"

TEST_FUNCTION_START(acb_modular_theta, state)
{
    slong iter;

    /* Test consistency with/without transform */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t t1, t2, t3, t4, t1b, t2b, t3b, t4b, z, tau;
        slong prec0, prec1, prec2, e0;

        acb_init(t1); acb_init(t2); acb_init(t3); acb_init(t4);
        acb_init(t1b); acb_init(t2b); acb_init(t3b); acb_init(t4b);
        acb_init(z);
        acb_init(tau);

        prec0 = 2 + n_randint(state, 2000);
        prec1 = 2 + n_randint(state, 2000);
        prec2 = 2 + n_randint(state, 2000);
        e0 = 1 + n_randint(state, 100);

        acb_randtest(tau, state, prec0, e0);
        acb_randtest(z, state, prec0, e0);
        arb_abs(acb_imagref(tau), acb_imagref(tau));

        if (n_randint(state, 3) == 0)
            arb_set_si(acb_realref(tau), -10 + n_randint(state, 20));
        if (n_randint(state, 3) == 0)
            arb_zero(acb_imagref(z));
        if (n_randint(state, 3) == 0)
            arb_zero(acb_realref(z));

        acb_modular_theta(t1, t2, t3, t4, z, tau, prec1);
        acb_modular_theta_notransform(t1b, t2b, t3b, t4b, z, tau, prec2);

        if (!acb_overlaps(t1, t1b) || !acb_overlaps(t2, t2b) ||
            !acb_overlaps(t3, t3b) || !acb_overlaps(t4, t4b))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("z = "); acb_printd(z, 25); flint_printf("\n\n");
            flint_printf("tau = "); acb_printd(tau, 25); flint_printf("\n\n");

            flint_printf("t1  = "); acb_printd(t1, 15); flint_printf("\n\n");
            flint_printf("t1b = "); acb_printd(t1b, 15); flint_printf("\n\n");

            flint_printf("t2  = "); acb_printd(t2, 15); flint_printf("\n\n");
            flint_printf("t2b = "); acb_printd(t2b, 15); flint_printf("\n\n");

            flint_printf("t3  = "); acb_printd(t3, 15); flint_printf("\n\n");
            flint_printf("t3b = "); acb_printd(t3b, 15); flint_printf("\n\n");

            flint_printf("t4  = "); acb_printd(t4, 15); flint_printf("\n\n");
            flint_printf("t4b = "); acb_printd(t4b, 15); flint_printf("\n\n");

            flint_abort();
        }

        acb_clear(t1); acb_clear(t2); acb_clear(t3); acb_clear(t4);
        acb_clear(t1b); acb_clear(t2b); acb_clear(t3b); acb_clear(t4b);
        acb_clear(z);
        acb_clear(tau);
    }

    TEST_FUNCTION_END(state);
}
