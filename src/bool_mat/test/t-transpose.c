/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "bool_mat.h"

TEST_FUNCTION_START(bool_mat_transpose, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n;
        bool_mat_t a, b;

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        bool_mat_init(a, m, n);
        bool_mat_init(b, n, m);

        bool_mat_randtest(a, state);
        bool_mat_randtest(b, state);

        bool_mat_transpose(b, a);

        /* involution */
        {
            bool_mat_t c;
            bool_mat_init(c, m, n);
            bool_mat_randtest(c, state);
            bool_mat_transpose(c, b);
            if (!bool_mat_equal(c, a))
            {
                flint_printf("FAIL (involution)\n");
                flint_printf("m = %wd, n = %wd\n", m, n);
                flint_abort();
            }
            bool_mat_clear(c);
        }

        /* aliasing */
        if (bool_mat_is_square(a))
        {
            bool_mat_transpose(a, a);

            if (!bool_mat_equal(a, b))
            {
                flint_printf("FAIL (aliasing)\n");
                flint_abort();
            }
        }

        bool_mat_clear(a);
        bool_mat_clear(b);
    }

    TEST_FUNCTION_END(state);
}
