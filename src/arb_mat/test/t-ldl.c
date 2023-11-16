/*
    Copyright (C) 2016 Arb authors

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq.h"
#include "fmpq_mat.h"
#include "arb_mat.h"

/* Defined in t-cho.c, t-inv_cho_precomp.c, t-inv_ldl_precomp.c, t-ldl.c,
   t-solve_cho_precomp.c, t-solve_ldl_precomp.c, t-spd_inv.c, t-spd_solve.c */
#ifndef _fmpq_mat_randtest_positive_semidefinite
#define _fmpq_mat_randtest_positive_semidefinite _fmpq_mat_randtest_positive_semidefinite
void
_fmpq_mat_randtest_positive_semidefinite(fmpq_mat_t mat, flint_rand_t state, flint_bitcnt_t bits)
{
    slong n;
    fmpq_mat_t R, RT;
    if (!fmpq_mat_is_square(mat)) flint_abort(); /* assert */
    n = fmpq_mat_nrows(mat);
    fmpq_mat_init(R, n, n);
    fmpq_mat_init(RT, n, n);
    fmpq_mat_randtest(R, state, bits);
    fmpq_mat_transpose(RT, R);
    fmpq_mat_mul(mat, R, RT);
    fmpq_mat_clear(R);
    fmpq_mat_clear(RT);
}
#endif

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

TEST_FUNCTION_START(arb_mat_ldl, state)
{
    slong iter;

    /* check special matrices */
    {
        slong n;
        for (n = 1; n < 10; n++)
        {
            slong lprec;
            arb_mat_t L, A;
            arb_mat_init(L, n, n);
            arb_mat_init(A, n, n);
            for (lprec = 2; lprec < 10; lprec++)
            {
                int result;
                slong prec;
                prec = 1 << lprec;

                /* zero */
                arb_mat_zero(A);
                result = arb_mat_ldl(L, A, prec);
                if (result)
                {
                    flint_printf("FAIL (zero):\n");
                    flint_printf("n = %wd, prec = %wd\n", n, prec);
                    flint_printf("L = \n"); arb_mat_printd(L, 15);
                    flint_printf("\n\n");
                }

                /* negative identity */
                arb_mat_one(A);
                arb_mat_neg(A, A);
                result = arb_mat_ldl(L, A, prec);
                if (result)
                {
                    flint_printf("FAIL (negative identity):\n");
                    flint_printf("n = %wd, prec = %wd\n", n, prec);
                    flint_printf("L = \n"); arb_mat_printd(L, 15);
                    flint_printf("\n\n");
                }

                /* identity */
                arb_mat_one(A);
                result = arb_mat_ldl(L, A, prec);
                if (!result || !arb_mat_equal(L, A))
                {
                    flint_printf("FAIL (identity):\n");
                    flint_printf("n = %wd, prec = %wd\n", n, prec);
                    flint_printf("L = \n"); arb_mat_printd(L, 15);
                    flint_printf("\n\n");
                }
            }
            arb_mat_clear(L);
            arb_mat_clear(A);
        }
    }

    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_mat_t Q;
        arb_mat_t A, L, D, U, T;
        slong n, qbits, prec;
        int q_invertible, r_invertible;

        n = n_randint(state, 8);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 202);

        fmpq_mat_init(Q, n, n);
        arb_mat_init(A, n, n);
        arb_mat_init(L, n, n);
        arb_mat_init(D, n, n);
        arb_mat_init(U, n, n);
        arb_mat_init(T, n, n);

        _fmpq_mat_randtest_positive_semidefinite(Q, state, qbits);
        q_invertible = fmpq_mat_is_invertible(Q);

        if (!q_invertible)
        {
            arb_mat_set_fmpq_mat(A, Q, prec);
            r_invertible = arb_mat_ldl(L, A, prec);
            if (r_invertible)
            {
                flint_printf("FAIL: matrix is singular over Q but not over R\n");
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("L = \n"); arb_mat_printd(L, 15); flint_printf("\n\n");
            }
        }
        else
        {
            /* now this must converge */
            while (1)
            {
                arb_mat_set_fmpq_mat(A, Q, prec);
                r_invertible = arb_mat_ldl(L, A, prec);
                if (r_invertible)
                {
                    break;
                }
                else
                {
                    if (prec > 10000)
                    {
                        flint_printf("FAIL: failed to converge at 10000 bits\n");
                        flint_printf("n = %wd, prec = %wd\n", n, prec);
                        flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                        flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                        flint_abort();
                    }
                    prec *= 2;
                }
            }

            /* multiply out the decomposition */
            {
                slong i;
                arb_mat_zero(D);
                arb_mat_transpose(U, L);
                for (i = 0; i < n; i++)
                {
                    arb_set(arb_mat_entry(D, i, i), arb_mat_entry(L, i, i));
                    arb_one(arb_mat_entry(L, i, i));
                    arb_one(arb_mat_entry(U, i, i));
                }
                arb_mat_mul(T, L, D, prec);
                arb_mat_mul(T, T, U, prec);
            }

            if (!arb_mat_contains_fmpq_mat(T, Q))
            {
                flint_printf("FAIL (containment, iter = %wd)\n", iter);
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("L = \n"); arb_mat_printd(L, 15); flint_printf("\n\n");
                flint_printf("U = \n"); arb_mat_printd(U, 15); flint_printf("\n\n");
                flint_printf("L*U = \n"); arb_mat_printd(T, 15); flint_printf("\n\n");

                flint_abort();
            }
        }

        fmpq_mat_clear(Q);
        arb_mat_clear(A);
        arb_mat_clear(L);
        arb_mat_clear(D);
        arb_mat_clear(U);
        arb_mat_clear(T);
    }

    TEST_FUNCTION_END(state);
}
