/*
    Copyright (C) 2012 Fredrik Johansson

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

void
_fmpq_mul_si_frac(fmpq_t res, const fmpq_t x, slong p, slong q)
{
    fmpq_t r;
    fmpq_init(r);
    fmpq_set_si(r, p, q);
    fmpq_mul(res, x, r);
    fmpq_clear(r);
}

void
_fmpq_max(fmpq_t p, const fmpq_t q, const fmpq_t r)
{
    if (fmpq_cmp(q, r) < 0)
        fmpq_set(p, r);
    else
        fmpq_set(p, q);
}

void
_fmpq_mat_scalar_mul_fmpq(fmpq_mat_t res,
        const fmpq_mat_t mat, const fmpq_t x)
{
    slong i, j;
    for (i = 0; i < mat->r; i++)
        for (j = 0; j < mat->c; j++)
            fmpq_mul(fmpq_mat_entry(res, i, j),
                     fmpq_mat_entry(mat, i, j), x);
}

void
_fmpq_mat_inf_norm(fmpq_t res, const fmpq_mat_t mat)
{
    fmpq_t s, q;
    slong i, j;

    fmpq_init(s);
    fmpq_init(q);

    fmpq_zero(res);
    for (i = 0; i < mat->r; i++)
    {
        fmpq_zero(s);
        for (j = 0; j < mat->c; j++)
        {
            fmpq_abs(q, fmpq_mat_entry(mat, i, j));
            fmpq_add(s, s, q);
        }
        _fmpq_max(res, res, s);
    }

    fmpq_clear(s);
    fmpq_clear(q);
}

TEST_FUNCTION_START(arb_mat_solve_preapprox, state)
{
    slong iter;

    for (iter = 0; iter < 2000 * 0.1 * flint_test_multiplier(); iter++)
    {
        fmpq_mat_t QA, QB, QR, QT, QE, QX;
        arb_mat_t A, B, R, T, X;
        slong n, m, qbits, prec;
        int q_invertible, r_approximable, r_approximable2;

        n = n_randint(state, 8);
        m = n_randint(state, 8);
        qbits = 1 + n_randint(state, 30);
        prec = 2 + n_randint(state, 200);

        fmpq_mat_init(QA, n, n);
        fmpq_mat_init(QB, n, m);
        fmpq_mat_init(QR, n, n);
        fmpq_mat_init(QT, n, m);
        fmpq_mat_init(QE, n, n);
        fmpq_mat_init(QX, n, m);

        arb_mat_init(A, n, n);
        arb_mat_init(B, n, m);
        arb_mat_init(R, n, n);
        arb_mat_init(T, n, m);
        arb_mat_init(X, n, m);

        /* Sample a square matrix QE with inf-norm less than 1. */
        {
            fmpq_t m;
            fmpq_init(m);
            fmpq_mat_randtest(QE, state, qbits);
            _fmpq_mat_inf_norm(m, QE);
            if (!fmpq_is_zero(m))
            {
                slong p, q;
                q = 64;
                p = n_randint(state, q);
                fmpq_inv(m, m);
                _fmpq_mul_si_frac(m, m, p, q);
                _fmpq_mat_scalar_mul_fmpq(QE, QE, m);
            }
            fmpq_clear(m);
        }

        /* Sample an unrestricted square matrix QR. */
        fmpq_mat_randtest(QR, state, qbits);

        /* Construct QA := inv(QR)(I - QE), so that I - QR*QA = QE */
        fmpq_mat_one(QA);
        fmpq_mat_sub(QA, QA, QE);
        q_invertible = fmpq_mat_solve_fraction_free(QA, QR, QA);

        /* Sample unrestricted matrices QB and QT. */
        fmpq_mat_randtest(QB, state, qbits);
        fmpq_mat_randtest(QT, state, qbits);

        arb_mat_set_fmpq_mat(A, QA, prec);
        arb_mat_set_fmpq_mat(B, QB, prec);
        arb_mat_set_fmpq_mat(R, QR, prec);
        arb_mat_set_fmpq_mat(T, QT, prec);

        if (q_invertible)
        {
            /* Construct QX := inv(QA)QB */
            fmpq_mat_solve_fraction_free(QX, QA, QB);

            r_approximable = arb_mat_solve_preapprox(X, A, B, R, T, prec);
            if (r_approximable)
            {
                if (!arb_mat_contains_fmpq_mat(X, QX))
                {
                    flint_printf("FAIL (containment, iter = %wd)\n", iter);
                    flint_printf("n = %wd, prec = %wd\n", n, prec);
                    flint_printf("\n");

                    flint_printf("QA = \n"); fmpq_mat_print(QA); flint_printf("\n\n");
                    flint_printf("QB = \n"); fmpq_mat_print(QB); flint_printf("\n\n");
                    flint_printf("QR = \n"); fmpq_mat_print(QR); flint_printf("\n\n");
                    flint_printf("QT = \n"); fmpq_mat_print(QT); flint_printf("\n\n");
                    flint_printf("QE = \n"); fmpq_mat_print(QE); flint_printf("\n\n");
                    flint_printf("QX = \n"); fmpq_mat_print(QX); flint_printf("\n\n");

                    flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                    flint_printf("B = \n"); arb_mat_printd(B, 15); flint_printf("\n\n");
                    flint_printf("R = \n"); arb_mat_printd(R, 15); flint_printf("\n\n");
                    flint_printf("T = \n"); arb_mat_printd(T, 15); flint_printf("\n\n");
                    flint_printf("X = \n"); arb_mat_printd(X, 15); flint_printf("\n\n");

                    flint_abort();
                }
            }

            /* test aliasing */
            if (n_randint(state, 2))
            {
                r_approximable2 = arb_mat_solve_preapprox(T, A, B, R, T, prec);
                if ((r_approximable != r_approximable2) ||
                    (r_approximable && !arb_mat_equal(X, T)))
                {
                    flint_printf("FAIL (aliasing (X, T))\n");
                    flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                    flint_printf("B = \n"); arb_mat_printd(B, 15); flint_printf("\n\n");
                    flint_printf("T = \n"); arb_mat_printd(T, 15); flint_printf("\n\n");
                    flint_printf("X = \n"); arb_mat_printd(X, 15); flint_printf("\n\n");
                    flint_abort();
                }
            }
            else
            {
                r_approximable2 = arb_mat_solve_preapprox(B, A, B, R, T, prec);
                if ((r_approximable != r_approximable2) ||
                    (r_approximable && !arb_mat_equal(X, B)))
                {
                    flint_printf("FAIL (aliasing (X, B))\n");
                    flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n\n");
                    flint_printf("B = \n"); arb_mat_printd(B, 15); flint_printf("\n\n");
                    flint_printf("T = \n"); arb_mat_printd(T, 15); flint_printf("\n\n");
                    flint_printf("X = \n"); arb_mat_printd(X, 15); flint_printf("\n\n");
                    flint_abort();
                }
            }
        }

        /* test a simplified case where A = R = I */
        {
            arb_mat_one(A);
            arb_mat_one(R);
            fmpq_mat_randtest(QB, state, qbits);
            fmpq_mat_randtest(QT, state, qbits);
            arb_mat_set_fmpq_mat(B, QB, prec);
            arb_mat_set_fmpq_mat(T, QT, prec);
            arb_mat_zero(X);
            r_approximable = arb_mat_solve_preapprox(X, A, B, R, T, prec);
            if (!r_approximable || !arb_mat_contains(X, B))
            {
                flint_printf("FAIL (special values A = R = I))\n");
                flint_printf("B = \n"); arb_mat_printd(B, 15); flint_printf("\n\n");
                flint_printf("T = \n"); arb_mat_printd(T, 15); flint_printf("\n\n");
                flint_printf("X = \n"); arb_mat_printd(X, 15); flint_printf("\n\n");
            }
        }

        /* test a simplified case where A = R = I and T = B */
        {
            arb_mat_one(A);
            arb_mat_one(R);
            fmpq_mat_randtest(QB, state, qbits);
            arb_mat_set_fmpq_mat(B, QB, prec);
            arb_mat_set(T, B);
            arb_mat_zero(X);
            r_approximable = arb_mat_solve_preapprox(X, A, B, R, T, prec);
            if (!r_approximable || !arb_mat_contains(X, B))
            {
                flint_printf("FAIL (special values A = R = I, T = B))\n");
                flint_printf("B = \n"); arb_mat_printd(B, 15); flint_printf("\n\n");
                flint_printf("X = \n"); arb_mat_printd(X, 15); flint_printf("\n\n");
            }
        }

        fmpq_mat_clear(QA);
        fmpq_mat_clear(QB);
        fmpq_mat_clear(QR);
        fmpq_mat_clear(QT);
        fmpq_mat_clear(QE);
        fmpq_mat_clear(QX);

        arb_mat_clear(A);
        arb_mat_clear(B);
        arb_mat_clear(R);
        arb_mat_clear(T);
        arb_mat_clear(X);
    }

    TEST_FUNCTION_END(state);
}
