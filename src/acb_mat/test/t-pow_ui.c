/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"
#include "acb_mat.h"

TEST_FUNCTION_START(acb_mat_pow_ui, state)
{
    slong iter;

    for (iter = 0; iter < 200 * flint_test_multiplier(); iter++)
    {
        fmpq_mat_t Q, R, S;
        acb_mat_t A, B;
        slong n, qbits, prec;
        ulong e, i;

        n = n_randint(state, 6);
        e = n_randint(state, 12);
        qbits = 1 + n_randint(state, 20);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_mat_init(R, n, n);
        fmpq_mat_init(S, n, n);
        acb_mat_init(A, n, n);
        acb_mat_init(B, n, n);

        fmpq_mat_randtest(Q, state, qbits);
        acb_mat_set_fmpq_mat(A, Q, prec);
        acb_mat_randtest(B, state, 100, 10);

        fmpq_mat_one(R);
        for (i = 0; i < e; i++)
        {
            fmpq_mat_mul(S, R, Q);
            fmpq_mat_swap(R, S);
        }

        acb_mat_pow_ui(B, A, e, prec);

        if (!acb_mat_contains_fmpq_mat(B, R))
        {
            flint_printf("FAIL (containment)\n\n");
            flint_printf("n = %wd, e = %wu, prec = %wd\n\n", n, e, prec);
            flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
            flint_printf("R = \n"); fmpq_mat_print(R); flint_printf("\n\n");
            flint_printf("B = \n"); acb_mat_printd(B, 15); flint_printf("\n\n");
            flint_abort();
        }

        /* aliasing */
        acb_mat_pow_ui(A, A, e, prec);

        if (!acb_mat_equal(A, B))
        {
            flint_printf("FAIL (aliasing)\n\n");
            flint_printf("n = %wd, e = %wu, prec = %wd\n\n", n, e, prec);
            flint_abort();
        }

        fmpq_mat_clear(Q);
        fmpq_mat_clear(R);
        fmpq_mat_clear(S);
        acb_mat_clear(A);
        acb_mat_clear(B);
    }

    TEST_FUNCTION_END(state);
}
