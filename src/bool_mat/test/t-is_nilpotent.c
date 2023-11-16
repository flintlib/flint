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

TEST_FUNCTION_START(bool_mat_is_nilpotent, state)
{
    slong iter;

    /* empty matrix is not nilpotent(?) */
    {
        bool_mat_t A;
        bool_mat_init(A, 0, 0);
        if (bool_mat_is_nilpotent(A))
        {
            flint_printf("FAIL (empty)\n");
            flint_abort();
        }
        bool_mat_clear(A);
    }

    /* detect nilpotency of matrices that are constructed to be nilpotent */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m;
        bool_mat_t A;

        m = n_randint(state, 10) + 1;
        bool_mat_init(A, m, m);
        bool_mat_randtest_nilpotent(A, state);

        if (!bool_mat_is_nilpotent(A))
        {
            flint_printf("FAIL (nilpotent by construction)\n");
            flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n");
            flint_abort();
        }

        bool_mat_clear(A);
    }

    /* check nilpotency by computing a matrix power */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m;
        bool_mat_t A, B;

        m = n_randint(state, 10) + 1;

        bool_mat_init(A, m, m);
        bool_mat_init(B, m, m);

        bool_mat_randtest(A, state);
        bool_mat_pow_ui(B, A, m);

        if (bool_mat_is_nilpotent(A) != !bool_mat_any(B))
        {
            flint_printf("FAIL (A^m)\n");
            flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n");
            flint_printf("B = \n"); bool_mat_print(B); flint_printf("\n");
            flint_abort();
        }

        bool_mat_clear(A);
        bool_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
