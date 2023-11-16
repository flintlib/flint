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
#include "arb_mat.h"

TEST_FUNCTION_START(arb_mat_sqr, state)
{
    slong iter;

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong n, qbits1, rbits1, rbits2;
        fmpq_mat_t A, B;
        arb_mat_t a, b, c;

        qbits1 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);

        n = n_randint(state, 10);

        fmpq_mat_init(A, n, n);
        fmpq_mat_init(B, n, n);

        arb_mat_init(a, n, n);
        arb_mat_init(b, n, n);
        arb_mat_init(c, n, n);

        fmpq_mat_randtest(A, state, qbits1);
        fmpq_mat_mul(B, A, A);

        arb_mat_set_fmpq_mat(a, A, rbits1);
        arb_mat_sqr(b, a, rbits2);

        if (!arb_mat_contains_fmpq_mat(b, B))
        {
            flint_printf("FAIL\n\n");
            flint_printf("n = %wd, bits2 = %wd\n", n, rbits2);

            flint_printf("A = "); fmpq_mat_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_mat_print(B); flint_printf("\n\n");

            flint_printf("a = "); arb_mat_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); arb_mat_printd(b, 15); flint_printf("\n\n");

            flint_abort();
        }

        /* test aliasing */
        arb_mat_set(c, a);
        arb_mat_sqr(c, c, rbits2);
        if (!arb_mat_equal(c, b))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_abort();
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);

        arb_mat_clear(a);
        arb_mat_clear(b);
        arb_mat_clear(c);
    }

    TEST_FUNCTION_END(state);
}
