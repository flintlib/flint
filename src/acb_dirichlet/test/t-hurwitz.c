/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_hurwitz, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, a, b, z1, z2;
        slong prec1, prec2;

        acb_init(s);
        acb_init(a);
        acb_init(b);
        acb_init(z1);
        acb_init(z2);

        prec1 = 2 + n_randint(state, 200);
        prec2 = 2 + n_randint(state, 200);

        acb_randtest(s, state, 1 + n_randint(state, 500), 3);
        acb_randtest(a, state, 1 + n_randint(state, 500), 3);
        acb_randtest(b, state, 1 + n_randint(state, 500), 10);
        acb_randtest(z1, state, 1 + n_randint(state, 500), 10);
        acb_randtest(z2, state, 1 + n_randint(state, 500), 10);

        if (n_randint(state, 2))
            acb_set_si(a, n_randint(state, 100) - 50);
        if (n_randint(state, 2))
            acb_set_si(s, n_randint(state, 100) - 50);

        /* test zeta(s,a) = 2^(-s) (zeta(s,a/2) + zeta(s,(a+1)/2)) */
        acb_mul_2exp_si(b, a, -1);
        acb_dirichlet_hurwitz(z2, s, b, prec2);
        acb_add_ui(b, a, 1, prec2);
        acb_mul_2exp_si(b, b, -1);
        acb_dirichlet_hurwitz(b, s, b, prec2);
        acb_add(z2, z2, b, prec2);
        acb_neg(b, s);
        acb_set_ui(z1, 2);
        acb_pow(b, z1, b, prec2);
        acb_mul(z2, z2, b, prec2);

        acb_dirichlet_hurwitz(z1, s, a, prec1);

        if (!acb_overlaps(z1, z2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("s = "); acb_printd(s, 30); flint_printf("\n\n");
            flint_printf("a = "); acb_printd(a, 30); flint_printf("\n\n");
            flint_printf("z1 = "); acb_printd(z1, 30); flint_printf("\n\n");
            flint_printf("z2 = "); acb_printd(z2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(s);
        acb_clear(a);
        acb_clear(b);
        acb_clear(z1);
        acb_clear(z2);
    }

    TEST_FUNCTION_END(state);
}
