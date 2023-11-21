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

TEST_FUNCTION_START(bool_mat_is_transitive, state)
{
    slong iter;

    /* special matrices */
    {
        slong n;
        for (n = 0; n < 10; n++)
        {
            bool_mat_t A;
            bool_mat_init(A, n, n);

            /* identity matrices are transitive */
            bool_mat_one(A);
            if (!bool_mat_is_transitive(A))
            {
                flint_printf("FAIL (identity matrix)\n");
                flint_abort();
            }

            /* square zero matrices are transitive */
            bool_mat_zero(A);
            if (!bool_mat_is_transitive(A))
            {
                flint_printf("FAIL (zero matrix)\n");
                flint_abort();
            }

            bool_mat_clear(A);
        }
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong n;
        bool_mat_t A;

        n = n_randint(state, 10);

        bool_mat_init(A, n, n);

        /* all square diagonal matrices are transitive */
        bool_mat_randtest_diagonal(A, state);
        if (!bool_mat_is_transitive(A))
        {
            flint_printf("FAIL (diagonal)\n");
            flint_printf("A:\n"); bool_mat_print(A); flint_printf("\n");
            flint_abort();
        }

        bool_mat_clear(A);
    }

    TEST_FUNCTION_END(state);
}
