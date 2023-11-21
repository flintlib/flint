/*
    Copyright (C) 2015 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"
#include "acb_mat.h"

TEST_FUNCTION_START(acb_mat_sqr, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong n, qbits1, rbits1, rbits2;
        fmpq_mat_t A, B;
        acb_mat_t a, b, c;

        qbits1 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        n = n_randint(state, 10);

        fmpq_mat_init(A, n, n);
        fmpq_mat_init(B, n, n);

        acb_mat_init(a, n, n);
        acb_mat_init(b, n, n);
        acb_mat_init(c, n, n);

        fmpq_mat_randtest(A, state, qbits1);
        fmpq_mat_mul(B, A, A);

        acb_mat_set_fmpq_mat(a, A, rbits1);
        acb_mat_sqr(b, a, rbits2);

        if (!acb_mat_contains_fmpq_mat(b, B))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n = %wd, bits2 = %wd\n", n, rbits2);

            flint_printf("A = "); fmpq_mat_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_mat_print(B); flint_printf("\n\n");

            flint_printf("a = "); acb_mat_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_mat_printd(b, 15); flint_printf("\n\n");

            flint_abort();
        }

        /* test aliasing */
        acb_mat_set(c, a);
        acb_mat_sqr(c, c, rbits2);
        if (!acb_mat_equal(c, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);

        acb_mat_clear(a);
        acb_mat_clear(b);
        acb_mat_clear(c);
    }

    TEST_FUNCTION_END(state);
}
