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

    flint_printf("pfq_series_direct....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        acb_poly_struct a[4], b[4];
        acb_poly_t z, s1, s2;
        slong i, p, q, len1, len2, n1, n2, prec1, prec2;
        int regularized;

        p = n_randint(state, 3);
        q = n_randint(state, 3);

        len1 = n_randint(state, 8);
        len2 = n_randint(state, 8);

        prec1 = 2 + n_randint(state, 400);
        prec2 = 2 + n_randint(state, 400);

        if (n_randint(state, 2))
            n1 = -1;
        else
            n1 = n_randint(state, 50);

        if (n_randint(state, 2))
            n2 = -1;
        else
            n2 = n_randint(state, 50);

        regularized = n_randint(state, 2);

        acb_poly_init(z);
        acb_poly_init(s1);
        acb_poly_init(s2);
        for (i = 0; i < p; i++)
            acb_poly_init(a + i);
        for (i = 0; i < q; i++)
            acb_poly_init(b + i);

        acb_poly_randtest(z, state, 1 + n_randint(state, 10), 1 + n_randint(state, 500), 10);

        for (i = 0; i < p; i++)
            acb_poly_randtest(a + i, state, 1 + n_randint(state, 10), 1 + n_randint(state, 500), 3);
        for (i = 0; i < q; i++)
            acb_poly_randtest(b + i, state, 1 + n_randint(state, 10), 1 + n_randint(state, 500), 3);

        acb_hypgeom_pfq_series_direct(s1, a, p, b, q, z, regularized, n1, len1, prec1);
        acb_hypgeom_pfq_series_direct(s2, a, p, b, q, z, regularized, n2, len2, prec2);

        acb_poly_truncate(s1, FLINT_MIN(len1, len2));
        acb_poly_truncate(s2, FLINT_MIN(len1, len2));

        if (!acb_poly_overlaps(s1, s2))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("iter = %wd\n", iter);
            flint_printf("n1 = %wd, n2 = %wd    prec1 = %wd, prec2 = %wd\n\n", n1, n2, prec1, prec2);
            flint_printf("p = %wd, q = %wd\n\n", p, q);
            flint_printf("z = "); acb_poly_printd(z, 15); flint_printf("\n\n");

            for (i = 0; i < p; i++)
            {
                flint_printf("a[%wd] = ", i); acb_poly_printd(a + i, 15); flint_printf("\n\n");
            }

            for (i = 0; i < q; i++)
            {
                flint_printf("b[%wd] = ", i); acb_poly_printd(b + i, 15); flint_printf("\n\n");
            }

            flint_printf("s1 = "); acb_poly_printd(s1, 15); flint_printf("\n\n");
            flint_printf("s2 = "); acb_poly_printd(s2, 15); flint_printf("\n\n");
            acb_poly_sub(s1, s1, s2, prec1);
            flint_printf("diff = "); acb_poly_printd(s1, 15); flint_printf("\n\n");
            flint_abort();
        }

        acb_poly_clear(z);
        acb_poly_clear(s1);
        acb_poly_clear(s2);
        for (i = 0; i < p; i++)
            acb_poly_clear(a + i);
        for (i = 0; i < q; i++)
            acb_poly_clear(b + i);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

