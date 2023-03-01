/*
    Copyright (C) 2014 Fredrik Johansson

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

    flint_printf("pfq_sum_fme....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 5000 * arb_test_multiplier(); iter++)
    {
        acb_ptr a, b;
        acb_t z, s1, s2, t1, t2;
        slong i, p, q, n, prec1, prec2;

        p = n_randint(state, 4);
        q = n_randint(state, 4);
        n = n_randint(state, 100);
        prec1 = 2 + n_randint(state, 500);
        prec2 = 2 + n_randint(state, 500);

        acb_init(z);
        acb_init(s1);
        acb_init(s2);
        acb_init(t1);
        acb_init(t2);

        acb_randtest(z, state, 1 + n_randint(state, 500), 1 + n_randint(state, 5));
        acb_randtest(s1, state, 1 + n_randint(state, 500), 1 + n_randint(state, 5));
        acb_randtest(t1, state, 1 + n_randint(state, 500), 1 + n_randint(state, 5));
        acb_randtest(s2, state, 1 + n_randint(state, 500), 1 + n_randint(state, 5));
        acb_randtest(t2, state, 1 + n_randint(state, 500), 1 + n_randint(state, 5));

        a = _acb_vec_init(p);
        b = _acb_vec_init(q);

        for (i = 0; i < p; i++)
            acb_randtest(a + i, state, 1 + n_randint(state, 100), 1 + n_randint(state, 5));
        for (i = 0; i < q; i++)
            acb_randtest(b + i, state, 1 + n_randint(state, 100), 1 + n_randint(state, 5));

        acb_hypgeom_pfq_sum_forward(s1, t1, a, p, b, q, z, n, prec1);
        acb_hypgeom_pfq_sum_fme(s2, t2, a, p, b, q, z, n, prec2);

        if (!acb_overlaps(s1, s2) || !acb_overlaps(t1, t2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("z = "); acb_print(a); flint_printf("\n\n");
            flint_printf("s1 = "); acb_print(s1); flint_printf("\n\n");
            flint_printf("s2 = "); acb_print(s2); flint_printf("\n\n");
            flint_printf("t1 = "); acb_print(t1); flint_printf("\n\n");
            flint_printf("t2 = "); acb_print(t2); flint_printf("\n\n");
            flint_abort();
        }

        _acb_vec_clear(a, p);
        _acb_vec_clear(b, q);

        acb_clear(z);
        acb_clear(s1);
        acb_clear(s2);
        acb_clear(t1);
        acb_clear(t2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
