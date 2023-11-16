/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"

TEST_FUNCTION_START(arb_mat_transpose, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n;
        arb_mat_t a, b, c;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        arb_mat_init(a, m, n);
        arb_mat_init(b, n, m);
        arb_mat_init(c, m, n);

        arb_mat_randtest(a, state, 2 + n_randint(state, 100), 10);
        arb_mat_randtest(b, state, 2 + n_randint(state, 100), 10);
        arb_mat_randtest(c, state, 2 + n_randint(state, 100), 10);

        arb_mat_transpose(b, a);
        arb_mat_transpose(c, b);

        if (!arb_mat_equal(c, a))
        {
            flint_printf("FAIL\n\n");
            flint_printf("m = %wd, n = %wd\n", m, n);
            flint_abort();
        }

        if (arb_mat_nrows(a) == arb_mat_ncols(a))
        {
            arb_mat_transpose(c, a);
            arb_mat_transpose(a, a);

            if (!arb_mat_equal(a, c))
            {
                flint_printf("FAIL (aliasing)\n\n");
                flint_abort();
            }
        }

        arb_mat_clear(a);
        arb_mat_clear(b);
        arb_mat_clear(c);
    }

    TEST_FUNCTION_END(state);
}
