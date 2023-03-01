/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_elliptic.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("inv_p....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        acb_t tau, z, w, pw;
        slong prec;

        acb_init(tau);
        acb_init(z);
        acb_init(w);
        acb_init(pw);

        prec = 2 + n_randint(state, 400);

        acb_randtest(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        acb_randtest(w, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        acb_randtest(tau, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        if (arf_sgn(arb_midref(acb_imagref(tau))) < 0)
            acb_neg(tau, tau);

        acb_elliptic_inv_p(w, z, tau, prec);
        acb_elliptic_p(pw, w, tau, prec);

        if (!acb_contains(pw, z))
        {
            flint_printf("FAIL (containment)\n");
            flint_printf("tau = "); acb_printd(tau, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w = "); acb_printd(w, 30); flint_printf("\n\n");
            flint_printf("pw = "); acb_printd(pw, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(tau);
        acb_clear(z);
        acb_clear(w);
        acb_clear(pw);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

