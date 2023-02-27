/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_modular.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("delta....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        acb_t tau, z1, z2;
        slong e0, prec0, prec1, prec2;

        acb_init(tau);
        acb_init(z1);
        acb_init(z2);

        e0 = 1 + n_randint(state, 10);
        prec0 = 2 + n_randint(state, 1000);
        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);

        acb_randtest(tau, state, prec0, e0);
        acb_randtest(z1, state, prec0, e0);
        acb_randtest(z2, state, prec0, e0);

        /* compare with eta */
        acb_modular_delta(z1, tau, prec1);

        acb_modular_eta(z2, tau, prec2);
        acb_pow_ui(z2, z2, 24, prec2);

        if (!acb_overlaps(z1, z2))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("tau = "); acb_printd(tau, 15); flint_printf("\n\n");
            flint_printf("z1 = "); acb_printd(z1, 15); flint_printf("\n\n");
            flint_printf("z2 = "); acb_printd(z2, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_modular_delta(tau, tau, prec1);

        if (!acb_overlaps(z1, tau))
        {
            flint_printf("FAIL (aliasing)\n");
            flint_printf("tau = "); acb_printd(tau, 15); flint_printf("\n\n");
            flint_printf("z1 = "); acb_printd(z1, 15); flint_printf("\n\n");
            flint_printf("z2 = "); acb_printd(z2, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(tau);
        acb_clear(z1);
        acb_clear(z2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

