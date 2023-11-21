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

TEST_FUNCTION_START(bool_mat_complement, state)
{
    slong iter;

    /* zero matrix */
    {
        slong m, n;
        for (m = 1; m < 10; m++)
        {
            for (n = 1; n < 10; n++)
            {
                bool_mat_t zero, C;
                bool_mat_init(zero, m, n);
                bool_mat_init(C, m, n);
                bool_mat_zero(zero);
                bool_mat_complement(C, zero);
                if (bool_mat_any(zero) || bool_mat_all(zero) ||
                    !bool_mat_any(C) || !bool_mat_all(C))
                {
                    flint_printf("FAIL (zero matrix)\n");
                    flint_abort();
                }
                bool_mat_clear(zero);
                bool_mat_clear(C);
            }
        }
    }

    /* identity matrix */
    {
        slong n;
        for (n = 2; n < 10; n++)
        {
            bool_mat_t one, C;
            bool_mat_init(one, n, n);
            bool_mat_init(C, n, n);
            bool_mat_one(one);
            bool_mat_complement(C, one);
            if (!bool_mat_any(one) || bool_mat_all(one) ||
                !bool_mat_any(C) || bool_mat_all(C))
            {
                flint_printf("FAIL (identity matrix)\n");
                flint_abort();
            }
            bool_mat_clear(one);
            bool_mat_clear(C);
        }
    }

    /* random matrices */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n;
        bool_mat_t A, C;

        m = n_randint(state, 10) + 1;
        n = n_randint(state, 10) + 1;

        bool_mat_init(A, m, n);
        bool_mat_init(C, m, n);

        bool_mat_randtest(A, state);
        bool_mat_randtest(C, state);
        bool_mat_complement(C, A);

        if ((bool_mat_all(A) && bool_mat_any(C)) ||
            (bool_mat_all(C) && bool_mat_any(A)))
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        /* involution */
        {
            bool_mat_t CC;
            bool_mat_init(CC, m, n);
            bool_mat_randtest(CC, state);
            bool_mat_complement(CC, C);
            if (!bool_mat_equal(A, CC))
            {
                flint_printf("FAIL (involution)\n");
                flint_abort();
            }
            bool_mat_clear(CC);
        }

        /* aliasing */
        bool_mat_complement(A, A);
        if (!bool_mat_equal(A, C))
        {
            flint_printf("FAIL (aliasing)\n");
            flint_abort();
        }

        bool_mat_clear(A);
        bool_mat_clear(C);
    }

    TEST_FUNCTION_END(state);
}
