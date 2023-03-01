/*
    Copyright (C) 2012 Fredrik Johansson

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

    flint_printf("inv....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000 * arb_test_multiplier(); iter++)
    {
        fmpq_mat_t Q, Qinv;
        acb_mat_t A, Ainv;
        slong n, qbits, prec;
        int q_invertible, r_invertible, r_invertible2;

        n = n_randint(state, 8);
        qbits = 1 + n_randint(state, 30);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_mat_init(Qinv, n, n);

        acb_mat_init(A, n, n);
        acb_mat_init(Ainv, n, n);

        fmpq_mat_randtest(Q, state, qbits);
        q_invertible = fmpq_mat_inv(Qinv, Q);

        if (!q_invertible)
        {
            acb_mat_set_fmpq_mat(A, Q, prec);
            r_invertible = acb_mat_inv(Ainv, A, prec);
            if (r_invertible)
            {
                flint_printf("FAIL: matrix is singular over Q but not over R\n");
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("Ainv = \n"); acb_mat_printd(Ainv, 15); flint_printf("\n\n");
                flint_abort();
            }
        }
        else
        {
            /* now this must converge */
            while (1)
            {
                acb_mat_set_fmpq_mat(A, Q, prec);
                r_invertible = acb_mat_inv(Ainv, A, prec);

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
                        flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                        flint_abort();
                    }
                    prec *= 2;
                }
            }

            if (!acb_mat_contains_fmpq_mat(Ainv, Qinv))
            {
                flint_printf("FAIL (containment, iter = %wd)\n", iter);
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("Qinv = \n"); fmpq_mat_print(Qinv); flint_printf("\n\n");

                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("Ainv = \n"); acb_mat_printd(Ainv, 15); flint_printf("\n\n");

                flint_abort();
            }

            /* test aliasing */
            r_invertible2 = acb_mat_inv(A, A, prec);
            if (!acb_mat_equal(A, Ainv) || r_invertible != r_invertible2)
            {
                flint_printf("FAIL (aliasing)\n");
                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("Ainv = \n"); acb_mat_printd(Ainv, 15); flint_printf("\n\n");
                flint_abort();
            }
        }

        fmpq_mat_clear(Q);
        fmpq_mat_clear(Qinv);
        acb_mat_clear(A);
        acb_mat_clear(Ainv);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
