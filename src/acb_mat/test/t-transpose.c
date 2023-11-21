/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"

TEST_FUNCTION_START(acb_mat_transpose, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n;
        acb_mat_t a, b, c;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        acb_mat_init(a, m, n);
        acb_mat_init(b, n, m);
        acb_mat_init(c, m, n);

        acb_mat_randtest(a, state, 2 + n_randint(state, 100), 10);
        acb_mat_randtest(b, state, 2 + n_randint(state, 100), 10);
        acb_mat_randtest(c, state, 2 + n_randint(state, 100), 10);

        acb_mat_transpose(b, a);
        acb_mat_transpose(c, b);

        if (!acb_mat_equal(c, a))
        {
            flint_printf("FAIL\n\n");
            flint_printf("m = %wd, n = %wd\n", m, n);
            flint_abort();
        }

        if (acb_mat_nrows(a) == acb_mat_ncols(a))
        {
            acb_mat_transpose(c, a);
            acb_mat_transpose(a, a);

            if (!acb_mat_equal(a, c))
            {
                flint_printf("FAIL (aliasing)\n\n");
                flint_abort();
            }
        }

        acb_mat_clear(a);
        acb_mat_clear(b);
        acb_mat_clear(c);
    }

    TEST_FUNCTION_END(state);
}
