/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mat.h"
#include "acb_mat.h"

TEST_FUNCTION_START(acb_mat_eig_multiple, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_mat_t A, R;
        acb_ptr E, F;
        acb_t b;
        slong i, j, n, prec, count, count2;
        int result;

        n = n_randint(state, 8);
        prec = 2 + n_randint(state, 200);

        acb_init(b);
        acb_mat_init(A, n, n);
        acb_mat_init(R, n, n);
        E = _acb_vec_init(n);
        F = _acb_vec_init(n);

        if (n_randint(state, 10) != 0)
        {
            for (i = 0; i < n; i++)
                acb_randtest(E + i, state, prec, 2);
        }
        else
        {
            /* Randomly repeat eigenvalues. */
            for (i = 0; i < n; i++)
            {
                if (i == 0 || n_randint(state, 2))
                    acb_randtest(E + i, state, prec, 2);
                else
                    acb_set(E + i, E + n_randint(state, i));
            }
        }

        if (n_randint(state, 2))
        {
            for (i = 0; i < n; i++)
                acb_get_mid(E + i, E + i);
        }

        acb_mat_randtest_eig(A, state, E, prec);
        acb_mat_approx_eig_qr(F, NULL, R, A, NULL, 0, prec);

        /* Perturb F further. */
        if (n_randint(state, 10) == 0)
        {
            for (i = 0; i < n; i++)
            {
                acb_randtest(b, state, prec, 1);
                acb_mul_2exp_si(b, b, -n_randint(state, prec));
                acb_add(F + i, F + i, b, prec);
            }
        }

        /* Perturb R further. */
        if (n_randint(state, 10) == 0)
        {
            j = n_randint(state, n);

            for (i = 0; i < n; i++)
            {
                acb_randtest(b, state, prec, 1);
                acb_mul_2exp_si(b, b, -10 - n_randint(state, prec));
                acb_add(acb_mat_entry(R, i, j), acb_mat_entry(R, i, j), b, prec);
            }
        }

        if (n_randint(state, 2))
            result = acb_mat_eig_multiple_rump(F, A, E, R, prec);
        else
            result = acb_mat_eig_multiple(F, A, E, R, prec);

        if (result)
        {
            count = 0;
            count2 = 0;

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (j == 0 || !acb_equal(F + j, F + j - 1))
                        count += acb_contains(F + j, E + i);
                }

                for (j = 0; j < n; j++)
                {
                    if (j == 0 || !acb_equal(F + j, F + j - 1))
                        count2 += acb_overlaps(F + j, E + i);
                }
            }

            if (count != n || count2 != n)
            {
                flint_printf("FAIL: count\n\n");
                flint_printf("A = \n"); acb_mat_printd(A, 20); flint_printf("\n\n");
                flint_printf("R = \n"); acb_mat_printd(R, 20); flint_printf("\n\n");
                flint_printf("count = %wd, count2 = %wd\n\n", count, count2);
                flint_printf("E = \n");
                for (j = 0; j < n; j++)
                {
                    acb_printd(E + j, 20);
                    flint_printf("\n");
                }
                flint_printf("F = \n");
                for (j = 0; j < n; j++)
                {
                    acb_printd(F + j, 20);
                    flint_printf("\n");
                }
                flint_abort();
            }
        }

        acb_mat_clear(A);
        acb_mat_clear(R);
        _acb_vec_clear(E, n);
        _acb_vec_clear(F, n);
        acb_clear(b);
    }

    /* Test convergence for DFT matrices */
    for (iter = 0; iter < 50 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_mat_t A, R, QC;
        acb_ptr E;
        acb_t t;
        fmpq_mat_t Q, Qinv;
        slong i, n, c0, c1, c2, c3;
        slong prec;
        int algorithm, result;

        n = n_randint(state, 30);
        algorithm = n_randint(state, 2);
        acb_mat_init(A, n, n);
        acb_mat_init(R, n, n);
        E = _acb_vec_init(n);
        acb_init(t);
        acb_mat_init(QC, n, n);
        fmpq_mat_init(Q, n, n);
        fmpq_mat_init(Qinv, n, n);

        /* The current algorithm is not robust enough. */
#if 0
        do {
            fmpq_mat_randtest(Q, state, 2 + n_randint(state, 10));
        } while (!fmpq_mat_inv(Qinv, Q));
#else
        fmpq_mat_one(Q);
        fmpq_mat_one(Qinv);
#endif

        for (prec = 32; ; prec *= 2)
        {
            if (prec > 10000)
            {
                flint_printf("FAIL: unsuccessful, prec > 10000\n\n");
                flint_printf("algorithm = %d, iter %wd\n\n", algorithm, iter);
                flint_abort();
            }

            acb_mat_dft(A, 0, prec);

#if 0
            acb_mat_set_fmpq_mat(QC, Q, prec);
            acb_mat_mul(A, A, QC, prec);
            acb_mat_set_fmpq_mat(QC, Qinv, prec);
            acb_mat_mul(A, QC, A, prec);
#endif

            acb_mat_approx_eig_qr(E, NULL, R, A, NULL, 0, prec);

            if (algorithm == 0)
                result = acb_mat_eig_multiple_rump(E, A, E, R, prec);
            else
                result = acb_mat_eig_multiple(E, A, E, R, prec);

            /* Verify the known eigenvalues + multiplicities */
            if (result)
            {
                c0 = c1 = c2 = c3 = 0;

                for (i = 0; i < n; i++)
                {
                    acb_set_d_d(t, 1.0, 0.0);
                    c0 += acb_contains(E + i, t);
                    acb_set_d_d(t, -1.0, 0.0);
                    c1 += acb_contains(E + i, t);
                    acb_set_d_d(t, 0.0, 1.0);
                    c2 += acb_contains(E + i, t);
                    acb_set_d_d(t, 0.0, -1.0);
                    c3 += acb_contains(E + i, t);
                }

                result = (n == 0 || (c0 == (n+4)/4 && c1 == (n+2)/4 && c2 == (n-1)/4 && c3 == (n+1)/4));
            }

            if (result)
                break;
        }

        acb_mat_clear(A);
        acb_mat_clear(R);
        acb_mat_clear(QC);
        _acb_vec_clear(E, n);
        acb_clear(t);
        fmpq_mat_clear(Q);
        fmpq_mat_clear(Qinv);
    }

    TEST_FUNCTION_END(state);
}
