/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_l_euler_product, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, t, u;
        dirichlet_group_t G;
        dirichlet_char_t chi;
        ulong q, k;
        slong prec;

        acb_init(s);
        acb_init(t);
        acb_init(u);

        q = 1 + n_randint(state, 50);
        prec = 2 + n_randint(state, 500);
        k = n_randint(state, n_euler_phi(q));

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);
        dirichlet_char_index(chi, G, k);

        if (n_randint(state, 2))
        {
            acb_set_ui(s, 2 + n_randint(state, prec + 50));
        }
        else
        {
            acb_randtest(s, state, 2 + n_randint(state, 200), 2);
            arb_abs(acb_realref(s), acb_realref(s));
            if (n_randint(state, 2))
                acb_add_ui(s, s, n_randtest(state), prec);
        }

        if (n_randint(state, 2))
            acb_dirichlet_l_hurwitz(t, s, NULL, G, chi, prec);
        else
            acb_dirichlet_l_euler_product(t, s, G, chi, prec * 1.5);

        acb_dirichlet_l_euler_product(u, s, G, chi, prec);

        if (!acb_overlaps(t, u))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter = %wd  q = %wu  k = %wu  prec = %wd\n\n", iter, q, k, prec);
            flint_printf("s = "); acb_printn(s, 100, 0); flint_printf("\n\n");
            flint_printf("t = "); acb_printn(t, 100, 0); flint_printf("\n\n");
            flint_printf("u = "); acb_printn(u, 100, 0); flint_printf("\n\n");
            flint_abort();
        }

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
        acb_clear(s);
        acb_clear(t);
        acb_clear(u);
    }

    TEST_FUNCTION_END(state);
}
