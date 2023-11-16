/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_elliptic.h"

TEST_FUNCTION_START(acb_elliptic_p_p_prime, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_struct pj[2];
        acb_t tau, z, p, pp, g2, g3, t;
        slong prec;

        acb_init(tau);
        acb_init(z);
        acb_init(p);
        acb_init(pp);
        acb_init(pj + 0);
        acb_init(pj + 1);
        acb_init(g2);
        acb_init(g3);
        acb_init(t);

        prec = 2 + n_randint(state, 400);

        acb_randtest(z, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        acb_randtest(tau, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        if (arf_sgn(arb_midref(acb_imagref(tau))) < 0)
            acb_neg(tau, tau);

        acb_elliptic_p(p, z, tau, prec);
        acb_elliptic_p_prime(pp, z, tau, prec);
        acb_elliptic_p_jet(pj, z, tau, 2, prec);

        if (!acb_overlaps(p, pj + 0) || !acb_overlaps(pp, pj + 1))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("tau = "); acb_printd(tau, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("p = "); acb_printd(p, 30); flint_printf("\n\n");
            flint_printf("pp = "); acb_printd(pp, 30); flint_printf("\n\n");
            flint_printf("pj0 = "); acb_printd(pj + 0, 30); flint_printf("\n\n");
            flint_printf("pj1 = "); acb_printd(pj + 1, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_elliptic_invariants(g2, g3, tau, prec);
        acb_pow_ui(pj + 0, pp, 2, prec);

        acb_mul(t, p, g2, prec);
        acb_add(t, t, g3, prec);
        acb_pow_ui(pj + 1, p, 3, prec);
        acb_mul_ui(pj + 1, pj + 1, 4, prec);
        acb_sub(pj + 1, pj + 1, t, prec);

        if (!acb_overlaps(pj + 0, pj + 1))
        {
            flint_printf("FAIL (check pp^2 = 4p^3-g2*p-g3)\n");
            flint_printf("tau = "); acb_printd(tau, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("p = "); acb_printd(p, 30); flint_printf("\n\n");
            flint_printf("pp = "); acb_printd(pp, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(tau);
        acb_clear(z);
        acb_clear(p);
        acb_clear(pp);
        acb_clear(pj + 0);
        acb_clear(pj + 1);
        acb_clear(g2);
        acb_clear(g3);
        acb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
