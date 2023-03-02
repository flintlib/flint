/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "acb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("bernoulli_poly_ui....");
    fflush(stdout);

    flint_randinit(state);

    /* test multiplication theorem */
    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        acb_t x, t, res1, res2;
        ulong n, m, k;
        slong prec;

        n = n_randint(state, 50);
        m = 1 + n_randint(state, 5);
        prec = 2 + n_randint(state, 200);

        acb_init(x);
        acb_init(t);
        acb_init(res1);
        acb_init(res2);

        acb_randtest(x, state, 2 + n_randint(state, 200), 20);
        acb_randtest(res1, state, 2 + n_randint(state, 200), 20);

        acb_mul_ui(t, x, m, prec);
        acb_bernoulli_poly_ui(res1, n, t, prec);

        acb_zero(res2);
        for (k = 0; k < m; k++)
        {
            acb_set_ui(t, k);
            acb_div_ui(t, t, m, prec);
            acb_add(t, t, x, prec);
            acb_bernoulli_poly_ui(t, n, t, prec);
            acb_add(res2, res2, t, prec);
        }

        if (n > 0)
        {
            arb_ui_pow_ui(acb_realref(t), m, n - 1, prec);
            acb_mul_arb(res2, res2, acb_realref(t), prec);
        }
        else
        {
            acb_div_ui(res2, res2, m, prec);
        }

        if (!acb_overlaps(res1, res2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("n = %wu, m = %wu\n\n", n, m);
            flint_printf("x = "); acb_printd(x, 15); flint_printf("\n\n");
            flint_printf("res1 = "); acb_printd(res1, 15); flint_printf("\n\n");
            flint_printf("res2 = "); acb_printd(res2, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(x);
        acb_clear(t);
        acb_clear(res1);
        acb_clear(res2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

