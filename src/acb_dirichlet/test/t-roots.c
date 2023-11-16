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

TEST_FUNCTION_START(acb_dirichlet_roots, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_t x, y, z;
        acb_dirichlet_roots_t roots;
        ulong n, k;
        slong prec;
        slong iter2;

        n = 1 + n_randint(state, 1000);
        prec = 2 + n_randint(state, 200);

        acb_init(x);
        acb_init(y);
        acb_init(z);

        acb_dirichlet_roots_init(roots, n, n_randtest(state), prec);
        acb_unit_root(y, n, prec);

        for (iter2 = 0; iter2 <= FLINT_MIN(n, 20); iter2++)
        {
            k = n_randint(state, 2 * n);
            acb_dirichlet_root(x, roots, k, prec);
            acb_pow_ui(z, y, k, prec);

            if (!acb_overlaps(x, z))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("iter = %wd  n = %wu  k = %wu  prec = %wd\n\n", iter, n, k, prec);
                flint_printf("x = "); acb_printn(x, 30, 0); flint_printf("\n\n");
                flint_printf("z = "); acb_printn(z, 30, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_dirichlet_roots_clear(roots);
        acb_clear(x);
        acb_clear(y);
        acb_clear(z);
    }

    TEST_FUNCTION_END(state);
}
