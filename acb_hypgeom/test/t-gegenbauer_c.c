/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("gegenbauer_c....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        acb_t n, m, z, t, u, res1, res2;
        slong prec1, prec2;

        acb_init(n);
        acb_init(m);
        acb_init(z);
        acb_init(t);
        acb_init(u);
        acb_init(res1);
        acb_init(res2);

        prec1 = 2 + n_randint(state, 300);
        prec2 = 2 + n_randint(state, 300);

        if (n_randint(state, 2))
        {
            acb_set_si(m, n_randint(state, 20) - 10);
            acb_set_si(n, n_randint(state, 20) - 10);
        }
        else
        {
            acb_randtest_param(n, state, 1 + n_randint(state, 400), 10);
            acb_randtest_param(m, state, 1 + n_randint(state, 400), 10);
        }

        acb_randtest_param(z, state, 1 + n_randint(state, 400), 10);

        acb_hypgeom_gegenbauer_c(res1, n, m, z, prec1);

        acb_one(t);
        acb_mul_2exp_si(t, t, -1);
        acb_sub(t, m, t, prec2);
        acb_hypgeom_jacobi_p(res2, n, t, t, z, prec2);
        acb_add_ui(t, t, 1, prec2);
        acb_rising(t, t, n, prec2);
        acb_div(res2, res2, t, prec2);
        acb_mul_2exp_si(t, m, 1);
        acb_rising(t, t, n, prec2);
        acb_mul(res2, res2, t, prec2);

        if (!acb_overlaps(res1, res2))
        {
            flint_printf("FAIL: consistency 1\n\n");
            flint_printf("iter = %wd, prec1 = %wd, prec2 = %wd\n\n", iter, prec1, prec2);
            flint_printf("n = "); acb_printd(n, 30); flint_printf("\n\n");
            flint_printf("m = "); acb_printd(m, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("res1 = "); acb_printd(res1, 30); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 30); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(n);
        acb_clear(m);
        acb_clear(z);
        acb_clear(t);
        acb_clear(u);
        acb_clear(res1);
        acb_clear(res2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

