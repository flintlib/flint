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
#include "acb_poly.h"

TEST_FUNCTION_START(acb_dirichlet_hurwitz_precomp, state)
{
    slong iter;

    for (iter = 0; iter < 500 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t s, a, z1, z2;
        ulong p, q;
        slong prec1, prec2, A, K, N, i;
        acb_dirichlet_hurwitz_precomp_t pre;
        int deflate;

        prec1 = 2 + n_randint(state, 100);
        prec2 = 2 + n_randint(state, 100);
        A = 1 + n_randint(state, 10);
        K = 1 + n_randint(state, 10);
        N = 1 + n_randint(state, 10);
        deflate = (n_randint(state, 3) == 0);

        acb_init(s);
        acb_init(a);
        acb_init(z1);
        acb_init(z2);

        acb_randtest(s, state, 1 + n_randint(state, 200), 2);

        acb_dirichlet_hurwitz_precomp_init(pre, s, deflate, A, K, N, prec1);

        for (i = 0; i < 10; i++)
        {
            q = 1 + n_randint(state, 1000);
            p = 1 + n_randint(state, q);

            acb_dirichlet_hurwitz_precomp_eval(z1, pre, p, q, prec1);

            acb_set_ui(a, p);
            acb_div_ui(a, a, q, prec2);
            if (deflate)
                _acb_poly_zeta_cpx_series(z2, s, a, 1, 1, prec2);
            else
                acb_hurwitz_zeta(z2, s, a, prec2);

            if (!acb_overlaps(z1, z2))
            {
                flint_printf("FAIL! (overlap)");
                flint_printf("s = "); acb_printn(s, 50, 0); flint_printf("\n\n");
                flint_printf("A = %wd  K = %wd  N = %wd\n\n", A, K, N);
                flint_printf("p = %wu  q = %wu\n\n", p, q);
                flint_printf("z1 = "); acb_printn(z1, 50, 0); flint_printf("\n\n");
                flint_printf("z2 = "); acb_printn(z2, 50, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_dirichlet_hurwitz_precomp_clear(pre);

        acb_clear(s);
        acb_clear(a);
        acb_clear(z1);
        acb_clear(z2);
    }

    TEST_FUNCTION_END(state);
}
