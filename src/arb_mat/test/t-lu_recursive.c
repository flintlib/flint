/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "perm.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "arb_mat.h"

/* Defined in t-cho.c, t-ldl.c, t-lu.c, t-lu_recursive.c */
#ifndef fmpq_mat_is_invertible
#define fmpq_mat_is_invertible fmpq_mat_is_invertible
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
#endif

TEST_FUNCTION_START(arb_mat_lu_recursive, state)
{
    slong iter;

    /* Dummy test with rectangular matrices. Rectangular matrices are
       not actually supported (the output may be bogus), but the algorithm
       should at least not crash. */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, prec;
        slong *perm;
        arb_mat_t A, LU;

        n = n_randint(state, 20);
        m = n_randint(state, 20);
        prec = 2 + n_randint(state, 200);

        arb_mat_init(A, n, m);
        arb_mat_init(LU, n, m);
        perm = _perm_init(n);

        arb_mat_randtest(A, state, prec, 10);

        if (n_randint(state, 2))
        {
            arb_mat_lu_recursive(perm, LU, A, prec);
        }
        else
        {
            arb_mat_set(LU, A);
            arb_mat_lu_recursive(perm, LU, LU, prec);
        }

        arb_mat_clear(A);
        arb_mat_clear(LU);
        _perm_clear(perm);
    }

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_mat_t Q;
        arb_mat_t A, LU, P, L, U, T;
        slong i, j, n, qbits, prec, *perm;
        int q_invertible, r_invertible;

        n = n_randint(state, 20);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 202);

        fmpq_mat_init(Q, n, n);
        arb_mat_init(A, n, n);
        arb_mat_init(LU, n, n);
        arb_mat_init(P, n, n);
        arb_mat_init(L, n, n);
        arb_mat_init(U, n, n);
        arb_mat_init(T, n, n);
        perm = _perm_init(n);

        fmpq_mat_randtest(Q, state, qbits);
        q_invertible = fmpq_mat_is_invertible(Q);

        if (!q_invertible)
        {
            arb_mat_set_fmpq_mat(A, Q, prec);
            r_invertible = arb_mat_lu_recursive(perm, LU, A, prec);
            if (r_invertible)
            {
                flint_printf("FAIL: matrix is singular over Q but not over R\n");
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("LU = \n"); arb_mat_printd(LU, 15); flint_printf("\n\n");
            }
        }
        else
        {
            /* now this must converge */
            while (1)
            {
                arb_mat_set_fmpq_mat(A, Q, prec);
                r_invertible = arb_mat_lu_recursive(perm, LU, A, prec);
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

            arb_mat_one(L);
            for (i = 0; i < n; i++)
                for (j = 0; j < i; j++)
                    arb_set(arb_mat_entry(L, i, j),
                        arb_mat_entry(LU, i, j));

            for (i = 0; i < n; i++)
                for (j = i; j < n; j++)
                    arb_set(arb_mat_entry(U, i, j),
                        arb_mat_entry(LU, i, j));

            for (i = 0; i < n; i++)
                arb_one(arb_mat_entry(P, perm[i], i));

            arb_mat_mul(T, P, L, prec);
            arb_mat_mul(T, T, U, prec);

            if (!arb_mat_contains_fmpq_mat(T, Q))
            {
                flint_printf("FAIL (containment, iter = %wd)\n", iter);
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("LU = \n"); arb_mat_printd(LU, 15); flint_printf("\n\n");
                flint_printf("L = \n"); arb_mat_printd(L, 15); flint_printf("\n\n");
                flint_printf("U = \n"); arb_mat_printd(U, 15); flint_printf("\n\n");
                flint_printf("P*L*U = \n"); arb_mat_printd(T, 15); flint_printf("\n\n");

                flint_abort();
            }
        }

        fmpq_mat_clear(Q);
        arb_mat_clear(A);
        arb_mat_clear(LU);
        arb_mat_clear(P);
        arb_mat_clear(L);
        arb_mat_clear(U);
        arb_mat_clear(T);
        _perm_clear(perm);
    }

    TEST_FUNCTION_END(state);
}
