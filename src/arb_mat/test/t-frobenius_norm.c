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

static void
_fmpq_mat_sum_of_squares(fmpq_t res, const fmpq_mat_t Q)
{
    slong i, j;
    fmpq_zero(res);
    for (i = 0; i < fmpq_mat_nrows(Q); i++)
    {
        for (j = 0; j < fmpq_mat_ncols(Q); j++)
        {
            fmpq_addmul(res, fmpq_mat_entry(Q, i, j), fmpq_mat_entry(Q, i, j));
        }
    }
}

TEST_FUNCTION_START(arb_mat_frobenius_norm, state)
{
    slong iter;

    /* compare to the exact rational norm */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_mat_t Q;
        fmpq_t q;
        arb_mat_t A;
        slong n, qbits, prec;

        n = n_randint(state, 8);
        qbits = 1 + n_randint(state, 100);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(Q, n, n);
        fmpq_init(q);

        arb_mat_init(A, n, n);

        fmpq_mat_randtest(Q, state, qbits);
        _fmpq_mat_sum_of_squares(q, Q);

        arb_mat_set_fmpq_mat(A, Q, prec);

        /* check that the arb interval contains the exact value */
        {
            arb_t a;
            arb_init(a);

            arb_mat_frobenius_norm(a, A, prec);
            arb_mul(a, a, a, prec);

            if (!arb_contains_fmpq(a, q))
            {
                flint_printf("FAIL (containment, iter = %wd)\n", iter);
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("frobenius_norm(Q)^2 = \n");
                fmpq_print(q); flint_printf("\n\n");

                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("frobenius_norm(A)^2 = \n");
                arb_printd(a, 15); flint_printf("\n\n");
                flint_printf("frobenius_norm(A)^2 = \n");
                arb_print(a); flint_printf("\n\n");

                flint_abort();
            }

            arb_clear(a);
        }

        /* check that the upper bound is not less than the exact value */
        {
            mag_t b;
            fmpq_t y;

            mag_init(b);
            fmpq_init(y);

            arb_mat_bound_frobenius_norm(b, A);
            mag_mul(b, b, b);
            mag_get_fmpq(y, b);

            if (fmpq_cmp(q, y) > 0)
            {
                flint_printf("FAIL (bound, iter = %wd)\n", iter);
                flint_printf("n = %wd, prec = %wd\n", n, prec);
                flint_printf("\n");

                flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n\n");
                flint_printf("frobenius_norm(Q)^2 = \n");
                fmpq_print(q); flint_printf("\n\n");

                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                flint_printf("bound_frobenius_norm(A)^2 = \n");
                mag_printd(b, 15); flint_printf("\n\n");
                flint_printf("bound_frobenius_norm(A)^2 = \n");
                mag_print(b); flint_printf("\n\n");

                flint_abort();
            }

            mag_clear(b);
            fmpq_clear(y);
        }

        fmpq_mat_clear(Q);
        fmpq_clear(q);
        arb_mat_clear(A);
    }

    /* check trace(A^T A) = frobenius_norm(A)^2 */
    for (iter = 0; iter < 10000 * 0.1 * flint_test_multiplier(); iter++)
    {
        slong m, n, prec;
        arb_mat_t A, AT, ATA;
        arb_t t;

        prec = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);

        arb_mat_init(A, m, n);
        arb_mat_init(AT, n, m);
        arb_mat_init(ATA, n, n);
        arb_init(t);

        arb_mat_randtest(A, state, 2 + n_randint(state, 100), 10);
        arb_mat_transpose(AT, A);
        arb_mat_mul(ATA, AT, A, prec);
        arb_mat_trace(t, ATA, prec);
        arb_sqrt(t, t, prec);

        /* check the norm bound */
        {
            mag_t low, frobenius;

            mag_init(low);
            arb_get_mag_lower(low, t);

            mag_init(frobenius);
            arb_mat_bound_frobenius_norm(frobenius, A);

            if (mag_cmp(low, frobenius) > 0)
            {
                flint_printf("FAIL (bound)\n", iter);
                flint_printf("m = %wd, n = %wd, prec = %wd\n", m, n, prec);
                flint_printf("\n");

                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");

                flint_printf("lower(sqrt(trace(A^T A))) = \n");
                mag_printd(low, 15); flint_printf("\n\n");

                flint_printf("bound_frobenius_norm(A) = \n");
                mag_printd(frobenius, 15); flint_printf("\n\n");

                flint_abort();
            }

            mag_clear(low);
            mag_clear(frobenius);
        }

        /* check the norm interval */
        {
            arb_t frobenius;

            arb_init(frobenius);
            arb_mat_frobenius_norm(frobenius, A, prec);

            if (!arb_overlaps(t, frobenius))
            {
                flint_printf("FAIL (overlap)\n", iter);
                flint_printf("m = %wd, n = %wd, prec = %wd\n", m, n, prec);
                flint_printf("\n");

                flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");

                flint_printf("sqrt(trace(A^T A)) = \n");
                arb_printd(t, 15); flint_printf("\n\n");

                flint_printf("frobenius_norm(A) = \n");
                arb_printd(frobenius, 15); flint_printf("\n\n");

                flint_abort();
            }

            arb_clear(frobenius);
        }

        arb_mat_clear(A);
        arb_mat_clear(AT);
        arb_mat_clear(ATA);
        arb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
