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

TEST_FUNCTION_START(bool_mat_is_diagonal, state)
{
    slong iter;

    /* all zero matrices are diagonal */
    {
        slong m, n;
        for (m = 0; m < 10; m++)
        {
            for (n = 0; n < 10; n++)
            {
                bool_mat_t zero;
                bool_mat_init(zero, m, n);
                bool_mat_zero(zero);
                if (!bool_mat_is_diagonal(zero))
                {
                    flint_printf("FAIL (zero matrix)\n");
                    flint_abort();
                }
                bool_mat_clear(zero);
            }
        }
    }

    /* all identity matrices are diagonal */
    {
        slong n;
        for (n = 0; n < 10; n++)
        {
            bool_mat_t one;
            bool_mat_init(one, n, n);
            bool_mat_one(one);
            if (!bool_mat_is_diagonal(one))
            {
                flint_printf("FAIL (identity matrix)\n");
                flint_abort();
            }
            bool_mat_clear(one);
        }
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n;
        bool_mat_t A;

        m = n_randint(state, 10) + 1;
        n = n_randint(state, 10) + 1;

        bool_mat_init(A, m, n);

        /* random diagonal */
        {
            bool_mat_randtest_diagonal(A, state);
            if (!bool_mat_is_diagonal(A))
            {
                flint_printf("FAIL (random diagonal)\n");
                flint_printf("A:\n"); bool_mat_print(A); flint_printf("\n");
                flint_abort();
            }
        }

        /* random non-diagonal */
        {
            slong i, j;
            bool_mat_randtest(A, state);
            i = n_randint(state, m);
            j = n_randint(state, n);
            if (i != j)
            {
                bool_mat_set_entry(A, i, j, 1);
                if (bool_mat_is_diagonal(A))
                {
                    flint_printf("FAIL (random non-diagonal)\n");
                    flint_printf("A:\n"); bool_mat_print(A); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        bool_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
