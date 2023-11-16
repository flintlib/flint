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

TEST_FUNCTION_START(acb_modular_j, state)
{
    slong iter;

    /* Test SL2Z invariance */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t tau1, tau2, z1, z2;
        slong e0, prec0, prec1, prec2;
        psl2z_t g;

        psl2z_init(g);
        acb_init(tau1);
        acb_init(tau2);
        acb_init(z1);
        acb_init(z2);

        e0 = 1 + n_randint(state, 100);
        prec0 = 2 + n_randint(state, 2000);
        prec1 = 2 + n_randint(state, 2000);
        prec2 = 2 + n_randint(state, 2000);

        acb_randtest(tau1, state, prec0, e0);
        acb_randtest(tau2, state, prec0, e0);
        acb_randtest(z1, state, prec0, e0);
        acb_randtest(z2, state, prec0, e0);

        psl2z_randtest(g, state, 1 + n_randint(state, 200));

        acb_modular_transform(tau2, g, tau1, prec0);

        acb_modular_j(z1, tau1, prec1);
        acb_modular_j(z2, tau2, prec2);

        if (!acb_overlaps(z1, z2))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("tau1 = "); acb_print(tau1); flint_printf("\n\n");
            flint_printf("tau2 = "); acb_print(tau2); flint_printf("\n\n");
            flint_printf("z1 = "); acb_print(z1); flint_printf("\n\n");
            flint_printf("z2 = "); acb_print(z2); flint_printf("\n\n");
            flint_abort();
        }

        acb_modular_j(tau1, tau1, prec2);

        if (!acb_overlaps(z1, tau1))
        {
            flint_printf("FAIL (aliasing)\n");
            flint_printf("tau1 = "); acb_print(tau1); flint_printf("\n\n");
            flint_printf("tau2 = "); acb_print(tau2); flint_printf("\n\n");
            flint_printf("z1 = "); acb_print(z1); flint_printf("\n\n");
            flint_printf("z2 = "); acb_print(z2); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(tau1);
        acb_clear(tau2);
        acb_clear(z1);
        acb_clear(z2);
        psl2z_clear(g);
    }

    /* Test special values */
    for (iter = 0; iter < 100 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t tau, z;
        slong prec;

        acb_init(tau);
        acb_init(z);

        prec = 2 + n_randint(state, 2000);

        acb_randtest(z, state, prec, 10);

        acb_onei(tau);
        acb_modular_j(z, tau, prec);
        acb_sub_ui(z, z, 1728, prec);

        if (!acb_contains_zero(z))
        {
            flint_printf("FAIL (value 1)\n");
            flint_printf("tau = "); acb_print(tau); flint_printf("\n\n");
            flint_printf("z = "); acb_print(z); flint_printf("\n\n");
            flint_abort();
        }

        acb_set_ui(tau, 2);
        acb_div_ui(tau, tau, 3, prec);
        acb_exp_pi_i(tau, tau, prec);

        acb_modular_j(z, tau, prec);

        if (!acb_contains_zero(z))
        {
            flint_printf("FAIL (value 2)\n");
            flint_printf("tau = "); acb_print(tau); flint_printf("\n\n");
            flint_printf("z = "); acb_print(z); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(tau);
        acb_clear(z);
    }

    TEST_FUNCTION_END(state);
}
