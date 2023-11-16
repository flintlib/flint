/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"
#include "arb_mat.h"

TEST_FUNCTION_START(arb_mat_solve, state)
{
    slong iter;

    for (iter = 0; iter < 100000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_mat_t Q, QX, QB;
        arb_mat_t A, X, B;
        slong n, m, qbits, prec;
        int q_invertible, r_invertible, r_invertible2;

        n = n_randint(state, 8);
        m = n_randint(state, 8);
        qbits = 1 + n_randint(state, 30);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_mat_init(QX, n, m);
        fmpq_mat_init(QB, n, m);

        arb_mat_init(A, n, n);
        arb_mat_init(X, n, m);
        arb_mat_init(B, n, m);

        fmpq_mat_randtest(Q, state, qbits);
        fmpq_mat_randtest(QB, state, qbits);

        q_invertible = fmpq_mat_solve_fraction_free(QX, Q, QB);

        if (!q_invertible)
        {
            arb_mat_set_fmpq_mat(A, Q, prec);
            r_invertible = arb_mat_solve(X, A, B, prec);
            if (r_invertible)
            {
                flint_printf("FAIL: matrix is singular over Q but not over R\n");
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("QX = \n"); fmpq_mat_print(QX); flint_printf("\n\n");
                flint_printf("QB = \n"); fmpq_mat_print(QB); flint_printf("\n\n");
                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                flint_abort();
            }
        }
        else
        {
            /* now this must converge */
            while (1)
            {
                arb_mat_set_fmpq_mat(A, Q, prec);
                arb_mat_set_fmpq_mat(B, QB, prec);

                r_invertible = arb_mat_solve(X, A, B, prec);
                if (r_invertible)
                {
                    break;
                }
                else
                {
                    if (prec > 10000)
                    {
                        flint_printf("FAIL: failed to converge at 10000 bits\n");
                        flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                        flint_printf("QX = \n"); fmpq_mat_print(QX); flint_printf("\n\n");
                        flint_printf("QB = \n"); fmpq_mat_print(QB); flint_printf("\n\n");
                        flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                        flint_abort();
                    }
                    prec *= 2;
                }
            }

            if (!arb_mat_contains_fmpq_mat(X, QX))
            {
                flint_printf("FAIL (containment, iter = %wd)\n", iter);
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("QB = \n"); fmpq_mat_print(QB); flint_printf("\n\n");
                flint_printf("QX = \n"); fmpq_mat_print(QX); flint_printf("\n\n");

                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("B = \n"); arb_mat_printd(B, 15); flint_printf("\n\n");
                flint_printf("X = \n"); arb_mat_printd(X, 15); flint_printf("\n\n");

                flint_abort();
            }

            /* test aliasing */
            r_invertible2 = arb_mat_solve(B, A, B, prec);
            if (!arb_mat_equal(X, B) || r_invertible != r_invertible2)
            {
                flint_printf("FAIL (aliasing)\n");
                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("B = \n"); arb_mat_printd(B, 15); flint_printf("\n\n");
                flint_printf("X = \n"); arb_mat_printd(X, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        fmpq_mat_clear(Q);
        fmpq_mat_clear(QB);
        fmpq_mat_clear(QX);
        arb_mat_clear(A);
        arb_mat_clear(B);
        arb_mat_clear(X);
    }

    TEST_FUNCTION_END(state);
}
