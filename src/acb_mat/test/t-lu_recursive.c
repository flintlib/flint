/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_mat.h"

int fmpq_mat_is_invertible(const fmpq_mat_t A)
{
    int r;
    fmpq_t t;
    fmpq_init(t);
    fmpq_mat_det(t, A);
    r = !fmpq_is_zero(t);
    fmpq_clear(t);
    return r;
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("lu_recursive....");
    fflush(stdout);

    flint_randinit(state);

    /* Dummy test with rectangular matrices. Rectangular matrices are
       not actually supported (the output may be bogus), but the algorithm
       should at least not crash. */
    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        slong m, n, prec;
        slong *perm;
        acb_mat_t A, LU;

        n = n_randint(state, 20);
        m = n_randint(state, 20);
        prec = 2 + n_randint(state, 200);

        acb_mat_init(A, n, m);
        acb_mat_init(LU, n, m);
        perm = _perm_init(n);

        acb_mat_randtest(A, state, prec, 10);

        if (n_randint(state, 2))
        {
            acb_mat_lu_recursive(perm, LU, A, prec);
        }
        else
        {
            acb_mat_set(LU, A);
            acb_mat_lu_recursive(perm, LU, LU, prec);
        }

        acb_mat_clear(A);
        acb_mat_clear(LU);
        _perm_clear(perm);
    }

    for (iter = 0; iter < 2000 * arb_test_multiplier(); iter++)
    {
        fmpq_mat_t Q;
        acb_mat_t A, LU, P, L, U, T;
        slong i, j, n, qbits, prec, *perm;
        int q_invertible, r_invertible;

        n = n_randint(state, 20);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 202);

        fmpq_mat_init(Q, n, n);
        acb_mat_init(A, n, n);
        acb_mat_init(LU, n, n);
        acb_mat_init(P, n, n);
        acb_mat_init(L, n, n);
        acb_mat_init(U, n, n);
        acb_mat_init(T, n, n);
        perm = _perm_init(n);

        fmpq_mat_randtest(Q, state, qbits);
        q_invertible = fmpq_mat_is_invertible(Q);

        if (!q_invertible)
        {
            acb_mat_set_fmpq_mat(A, Q, prec);
            r_invertible = acb_mat_lu_recursive(perm, LU, A, prec);
            if (r_invertible)
            {
                flint_printf("FAIL: matrix is singular over Q but not over R\n");
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("LU = \n"); acb_mat_printd(LU, 15); flint_printf("\n\n");
            }
        }
        else
        {
            /* now this must converge */
            while (1)
            {
                acb_mat_set_fmpq_mat(A, Q, prec);
                r_invertible = acb_mat_lu_recursive(perm, LU, A, prec);
                if (r_invertible)
                {
                    break;
                }
                else
                {
                    if (prec > 10000)
                    {
                        flint_printf("FAIL: failed to converge at 10000 bits\n");
                        flint_abort();
                    }
                    prec *= 2;
                }
            }

            acb_mat_one(L);
            for (i = 0; i < n; i++)
                for (j = 0; j < i; j++)
                    acb_set(acb_mat_entry(L, i, j),
                        acb_mat_entry(LU, i, j));

            for (i = 0; i < n; i++)
                for (j = i; j < n; j++)
                    acb_set(acb_mat_entry(U, i, j),
                        acb_mat_entry(LU, i, j));

            for (i = 0; i < n; i++)
                acb_one(acb_mat_entry(P, perm[i], i));

            acb_mat_mul(T, P, L, prec);
            acb_mat_mul(T, T, U, prec);

            if (!acb_mat_contains_fmpq_mat(T, Q))
            {
                flint_printf("FAIL (containment, iter = %wd)\n", iter);
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("A = \n"); acb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("LU = \n"); acb_mat_printd(LU, 15); flint_printf("\n\n");
                flint_printf("L = \n"); acb_mat_printd(L, 15); flint_printf("\n\n");
                flint_printf("U = \n"); acb_mat_printd(U, 15); flint_printf("\n\n");
                flint_printf("P*L*U = \n"); acb_mat_printd(T, 15); flint_printf("\n\n");

                flint_abort();
            }
        }

        fmpq_mat_clear(Q);
        acb_mat_clear(A);
        acb_mat_clear(LU);
        acb_mat_clear(P);
        acb_mat_clear(L);
        acb_mat_clear(U);
        acb_mat_clear(T);
        _perm_clear(perm);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
