/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"
#include "acb_mat.h"

TEST_FUNCTION_START(acb_mat_mul_threaded, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, k, qbits1, qbits2, rbits1, rbits2, rbits3;
        fmpq_mat_t A, B, C;
        acb_mat_t a, b, c, d;

        flint_set_num_threads(1 + n_randint(state, 5));

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        k = n_randint(state, 10);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, n, k);
        fmpq_mat_init(C, m, k);

        acb_mat_init(a, m, n);
        acb_mat_init(b, n, k);
        acb_mat_init(c, m, k);
        acb_mat_init(d, m, k);

        fmpq_mat_randtest(A, state, qbits1);
        fmpq_mat_randtest(B, state, qbits2);
        fmpq_mat_mul(C, A, B);

        acb_mat_set_fmpq_mat(a, A, rbits1);
        acb_mat_set_fmpq_mat(b, B, rbits2);
        acb_mat_mul_threaded(c, a, b, rbits3);

        if (!acb_mat_contains_fmpq_mat(c, C))
        {
            flint_printf("FAIL\n\n");
            flint_printf("threads = %d, m = %wd, n = %wd, k = %wd, bits3 = %wd\n",
                flint_get_num_threads(), m, n, k, rbits3);

            flint_printf("A = "); fmpq_mat_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_mat_print(B); flint_printf("\n\n");
            flint_printf("C = "); fmpq_mat_print(C); flint_printf("\n\n");

            flint_printf("a = "); acb_mat_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_mat_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_mat_printd(c, 15); flint_printf("\n\n");

            flint_abort();
        }

        /* test aliasing with a */
        if (acb_mat_nrows(a) == acb_mat_nrows(c) &&
            acb_mat_ncols(a) == acb_mat_ncols(c))
        {
            acb_mat_set(d, a);
            acb_mat_mul_threaded(d, d, b, rbits3);
            if (!acb_mat_equal(d, c))
            {
                flint_printf("FAIL (aliasing 1)\n\n");
                flint_abort();
            }
        }

        /* test aliasing with b */
        if (acb_mat_nrows(b) == acb_mat_nrows(c) &&
            acb_mat_ncols(b) == acb_mat_ncols(c))
        {
            acb_mat_set(d, b);
            acb_mat_mul_threaded(d, a, d, rbits3);
            if (!acb_mat_equal(d, c))
            {
                flint_printf("FAIL (aliasing 2)\n\n");
                flint_abort();
            }
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);

        acb_mat_clear(a);
        acb_mat_clear(b);
        acb_mat_clear(c);
        acb_mat_clear(d);
    }

    TEST_FUNCTION_END(state);
}
