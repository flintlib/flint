/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb.h"
#include "acb_elliptic.h"

TEST_FUNCTION_START(acb_elliptic_invariants, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t tau, e1, e2, e3, g2, g3, h2, h3;
        slong prec;

        acb_init(tau);
        acb_init(e1);
        acb_init(e2);
        acb_init(e3);
        acb_init(g2);
        acb_init(g3);
        acb_init(h2);
        acb_init(h3);

        prec = 2 + n_randint(state, 400);

        acb_randtest(tau, state, 1 + n_randint(state, 200), 1 + n_randint(state, 10));
        if (arf_sgn(arb_midref(acb_imagref(tau))) < 0)
            acb_neg(tau, tau);

        acb_elliptic_roots(e1, e2, e3, tau, prec);
        acb_elliptic_invariants(g2, g3, tau, prec);

        acb_mul(h2, e1, e1, prec);
        acb_addmul(h2, e2, e2, prec);
        acb_addmul(h2, e3, e3, prec);
        acb_mul_2exp_si(h2, h2, 1);

        acb_mul(h3, e1, e2, prec);
        acb_mul(h3, h3, e3, prec);
        acb_mul_2exp_si(h3, h3, 2);

        if (!acb_overlaps(g2, h2) || !acb_overlaps(g3, h3))
        {
            flint_printf("FAIL (overlap)\n");
            flint_printf("tau = "); acb_printd(tau, 30); flint_printf("\n\n");
            flint_printf("g2 = "); acb_printd(g2, 30); flint_printf("\n\n");
            flint_printf("g3 = "); acb_printd(g3, 30); flint_printf("\n\n");
            flint_printf("h2 = "); acb_printd(h2, 30); flint_printf("\n\n");
            flint_printf("h3 = "); acb_printd(h3, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(tau);
        acb_clear(e1);
        acb_clear(e2);
        acb_clear(e3);
        acb_clear(g2);
        acb_clear(g3);
        acb_clear(h2);
        acb_clear(h3);
    }

    TEST_FUNCTION_END(state);
}
