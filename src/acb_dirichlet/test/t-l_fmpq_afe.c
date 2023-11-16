/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "acb_dirichlet.h"

TEST_FUNCTION_START(acb_dirichlet_l_fmpq_afe, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_t x;
        acb_t s, t, u;
        dirichlet_group_t G;
        dirichlet_char_t chi;
        ulong q, k;
        slong prec;

        fmpq_init(x);
        acb_init(s);
        acb_init(t);
        acb_init(u);

        q = 1 + n_randint(state, 50);
        prec = 2 + n_randint(state, 400);
        k = n_randint(state, n_euler_phi(q));

        dirichlet_group_init(G, q);
        dirichlet_char_init(chi, G);
        dirichlet_char_index(chi, G, k);

        if (dirichlet_char_is_primitive(G, chi))
        {
            if (n_randint(state, 2))
                fmpq_set_si(x, -8 + (slong) n_randint(state, 8), 1);
            else
                fmpq_randtest(x, state, 2 + n_randint(state, 8));

            acb_set_fmpq(s, x, prec);

            if (n_randint(state, 2))
                acb_dirichlet_l_hurwitz(t, s, NULL, G, chi, prec);
            else
                acb_dirichlet_l(t, s, G, chi, prec);

            acb_dirichlet_l_fmpq_afe(u, x, G, chi, prec);

            if (!acb_overlaps(t, u))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("iter = %wd  q = %wu  k = %wu  prec = %wd\n\n", iter, q, k, prec);
                flint_printf("x = "); fmpq_print(x); flint_printf("\n\n");
                flint_printf("s = "); acb_printn(s, 100, 0); flint_printf("\n\n");
                flint_printf("t = "); acb_printn(t, 100, 0); flint_printf("\n\n");
                flint_printf("u = "); acb_printn(u, 100, 0); flint_printf("\n\n");
                acb_sub(t, t, u, prec);
                flint_printf("t - u = "); acb_printd(t, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        dirichlet_char_clear(chi);
        dirichlet_group_clear(G);
        fmpq_clear(x);
        acb_clear(s);
        acb_clear(t);
        acb_clear(u);
    }

    TEST_FUNCTION_END(state);
}
