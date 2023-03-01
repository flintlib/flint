/*
    Copyright (C) 2018 arbguest

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("solve_precond....");
    fflush(stdout);

    flint_randinit(state);

    /* test random matrices, to test complex solving */
    for (iter = 0; iter < 2000 * arb_test_multiplier(); iter++)
    {
        acb_mat_t A, X, B, Y;
        slong n, m, prec;

        n = n_randint(state, 10);
        m = n_randint(state, 10);

        acb_mat_init(A, n, n);
        acb_mat_init(X, n, m);
        acb_mat_init(B, n, m);
        acb_mat_init(Y, n, m);

        prec = 2 + n_randint(state, 200);

        acb_mat_randtest(A, state, 200, 1 + n_randint(state, 10));
        acb_mat_randtest(X, state, 200, 1 + n_randint(state, 10));
        acb_mat_randtest(Y, state, 200, 1 + n_randint(state, 10));

        acb_mat_mul(B, A, X, prec);

        if (acb_mat_solve_precond(Y, A, B, prec))
        {
            if (!acb_mat_contains(Y, X))
            {
                flint_printf("FAIL: containment\n");
                flint_printf("A = \n"); acb_mat_printd(A, 30); flint_printf("\n\n");
                flint_printf("X = \n"); acb_mat_printd(A, 30); flint_printf("\n\n");
                flint_printf("B = \n"); acb_mat_printd(A, 30); flint_printf("\n\n");
                flint_printf("Y = \n"); acb_mat_printd(A, 30); flint_printf("\n\n");
                flint_abort();
            }
        }

        acb_mat_clear(A);
        acb_mat_clear(X);
        acb_mat_clear(B);
        acb_mat_clear(Y);
    }

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        fmpq_mat_t Q, QX, QB;
        acb_mat_t A, X, B;
        slong n, m, qbits, prec;
        int q_invertible, r_invertible, r_invertible2;

        n = n_randint(state, 8);
        m = n_randint(state, 8);
        qbits = 1 + n_randint(state, 30);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_mat_init(QX, n, m);
        fmpq_mat_init(QB, n, m);

        acb_mat_init(A, n, n);
        acb_mat_init(X, n, m);
        acb_mat_init(B, n, m);

        fmpq_mat_randtest(Q, state, qbits);
        fmpq_mat_randtest(QB, state, qbits);

        q_invertible = fmpq_mat_solve_fraction_free(QX, Q, QB);

        if (!q_invertible)
        {
            acb_mat_set_fmpq_mat(A, Q, prec);
            r_invertible = acb_mat_solve_precond(X, A, B, prec);
            if (r_invertible)
            {
                flint_printf("FAIL: matrix is singular over Q but not over C\n");
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("QX = \n"); fmpq_mat_print(QX); flint_printf("\n\n");
                flint_printf("QB = \n"); fmpq_mat_print(QB); flint_printf("\n\n");
                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                flint_abort();
            }
        }
        else
        {
            /* now this must converge */
            while (1)
            {
                acb_mat_set_fmpq_mat(A, Q, prec);
                acb_mat_set_fmpq_mat(B, QB, prec);

                r_invertible = acb_mat_solve_precond(X, A, B, prec);
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
                        flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                        flint_abort();
                    }
                    prec *= 2;
                }
            }

            if (!acb_mat_contains_fmpq_mat(X, QX))
            {
                flint_printf("FAIL (containment, iter = %wd)\n", iter);
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("QB = \n"); fmpq_mat_print(QB); flint_printf("\n\n");
                flint_printf("QX = \n"); fmpq_mat_print(QX); flint_printf("\n\n");

                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("B = \n"); acb_mat_printd(B, 15); flint_printf("\n\n");
                flint_printf("X = \n"); acb_mat_printd(X, 15); flint_printf("\n\n");

                flint_abort();
            }

            /* test aliasing */
            r_invertible2 = acb_mat_solve_precond(B, A, B, prec);
            if (!acb_mat_equal(X, B) || r_invertible != r_invertible2)
            {
                flint_printf("FAIL (aliasing)\n");
                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("B = \n"); acb_mat_printd(B, 15); flint_printf("\n\n");
                flint_printf("X = \n"); acb_mat_printd(X, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        fmpq_mat_clear(Q);
        fmpq_mat_clear(QB);
        fmpq_mat_clear(QX);
        acb_mat_clear(A);
        acb_mat_clear(B);
        acb_mat_clear(X);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
