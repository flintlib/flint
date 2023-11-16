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

TEST_FUNCTION_START(bool_mat_trace, state)
{
    slong iter;

    /* trace(empty) == 0 */
    {
        bool_mat_t A;
        bool_mat_init(A, 0, 0);
        bool_mat_one(A);
        if (bool_mat_trace(A) != 0)
        {
            flint_printf("FAIL (empty)\n");
            flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n\n");
            flint_abort();
        }
        bool_mat_clear(A);
    }

    /* trace(zero) == 0 */
    {
        slong n;
        bool_mat_t A;
        for (n = 1; n < 10; n++)
        {
            bool_mat_init(A, n, n);
            bool_mat_zero(A);
            if (bool_mat_trace(A) != 0)
            {
                flint_printf("FAIL (zero)\n");
                flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n\n");
                flint_abort();
            }
            bool_mat_clear(A);
        }
    }

    /* trace(one) == 1 */
    {
        slong n;
        bool_mat_t A;
        for (n = 1; n < 10; n++)
        {
            bool_mat_init(A, n, n);
            bool_mat_one(A);
            if (bool_mat_trace(A) != 1)
            {
                flint_printf("FAIL (one)\n");
                flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n\n");
                flint_abort();
            }
            bool_mat_clear(A);
        }
    }

    /* traces of random matrices with modified diagonal */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong n, i;
        bool_mat_t A;

        n = n_randint(state, 10) + 1;
        bool_mat_init(A, n, n);
        bool_mat_randtest(A, state);

        i = (slong) n_randint(state, n);
        bool_mat_set_entry(A, i, i, 1);
        if (bool_mat_trace(A) != 1)
        {
            flint_printf("FAIL (diagonal has a non-zero entry)\n");
            flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n\n");
            flint_abort();
        }

        for (i = 0; i < n; i++)
        {
            bool_mat_set_entry(A, i, i, 0);
        }
        if (bool_mat_trace(A) != 0)
        {
            flint_printf("FAIL (diagonal is zero)\n");
            flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n\n");
            flint_abort();
        }

        bool_mat_clear(A);
    }

    /* trace(A + B) == trace(A) + trace(B) */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong n;
        bool_mat_t A, B, C;

        n = n_randint(state, 10);

        bool_mat_init(A, n, n);
        bool_mat_init(B, n, n);
        bool_mat_init(C, n, n);

        bool_mat_randtest(A, state);
        bool_mat_randtest(B, state);
        bool_mat_add(C, A, B);

        if (bool_mat_trace(C) != (bool_mat_trace(A) | bool_mat_trace(B)))
        {
            flint_printf("FAIL (trace(A+B) == trace(A) | trace(B))\n");
            flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n\n");
            flint_printf("B = \n"); bool_mat_print(B); flint_printf("\n\n");
            flint_printf("A+B = \n"); bool_mat_print(C); flint_printf("\n\n");
            flint_abort();
        }

        bool_mat_clear(A);
        bool_mat_clear(B);
        bool_mat_clear(C);
    }

    /* trace(A*B) == trace(B*A) */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong n, m;
        bool_mat_t A, B, AB, BA;

        n = n_randint(state, 10);
        m = n_randint(state, 10);

        bool_mat_init(A, n, m);
        bool_mat_init(B, m, n);
        bool_mat_init(AB, n, n);
        bool_mat_init(BA, m, m);

        bool_mat_randtest(A, state);
        bool_mat_randtest(B, state);
        bool_mat_mul(AB, A, B);
        bool_mat_mul(BA, B, A);

        if (bool_mat_trace(AB) != bool_mat_trace(BA))
        {
            flint_printf("FAIL (trace(AB) == trace(BA))\n");
            flint_printf("A = \n"); bool_mat_print(A); flint_printf("\n\n");
            flint_printf("B = \n"); bool_mat_print(B); flint_printf("\n\n");
            flint_printf("AB = \n"); bool_mat_print(AB); flint_printf("\n\n");
            flint_printf("BA = \n"); bool_mat_print(BA); flint_printf("\n\n");
            flint_abort();
        }

        bool_mat_clear(A);
        bool_mat_clear(B);
        bool_mat_clear(AB);
        bool_mat_clear(BA);
    }

    TEST_FUNCTION_END(state);
}
