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
#include "acb_poly.h"
#include "acb_mat.h"

TEST_FUNCTION_START(acb_mat_eig_simple, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_mat_t A, L, R, LAR, D;
        acb_ptr E, F;
        acb_t b;
        slong i, j, n, prec, count, count2;
        int result, algorithm;

        n = n_randint(state, 8);
        prec = 2 + n_randint(state, 200);
        algorithm = n_randint(state, 3);

        acb_mat_init(A, n, n);
        acb_mat_init(L, n, n);
        acb_mat_init(R, n, n);
        acb_mat_init(LAR, n, n);
        acb_mat_init(D, n, n);
        acb_init(b);
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
        {
            if (algorithm == 0)
                result = acb_mat_eig_simple(F, L, R, A, E, R, prec);
            else if (algorithm == 1)
                result = acb_mat_eig_simple_rump(F, L, R, A, E, R, prec);
            else
                result = acb_mat_eig_simple_vdhoeven_mourrain(F, L, R, A, E, R, prec);
        }
        else
        {
            int r1, r2;
            if (algorithm == 0)
            {
                r1 = acb_mat_eig_simple(F, L, NULL, A, E, R, prec);
                r2 = acb_mat_eig_simple(F, NULL, R, A, E, R, prec);
            }
            else if (algorithm == 1)
            {
                r1 = acb_mat_eig_simple_rump(F, L, NULL, A, E, R, prec);
                r2 = acb_mat_eig_simple_rump(F, NULL, R, A, E, R, prec);
            }
            else
            {
                r1 = acb_mat_eig_simple_vdhoeven_mourrain(F, L, NULL, A, E, R, prec);
                r2 = acb_mat_eig_simple_vdhoeven_mourrain(F, NULL, R, A, E, R, prec);
            }
            result = n_randint(state, 2) ? r1 : r2;
        }

        acb_mat_mul(LAR, L, A, prec);
        acb_mat_mul(LAR, LAR, R, prec);
        for (i = 0; i < n; i++)
            acb_set(acb_mat_entry(D, i, i), F + i);

        if (!acb_mat_overlaps(LAR, D))
        {
            flint_printf("FAIL: overlap\n\n");
            flint_printf("algorithm = %d\n\n", algorithm);
            flint_printf("A = \n"); acb_mat_printd(A, 20); flint_printf("\n\n");
            flint_printf("R = \n"); acb_mat_printd(R, 20); flint_printf("\n\n");
            flint_printf("L = \n"); acb_mat_printd(L, 20); flint_printf("\n\n");
            flint_printf("D = \n"); acb_mat_printd(D, 20); flint_printf("\n\n");
            flint_printf("LAR = \n"); acb_mat_printd(LAR, 20); flint_printf("\n\n");
            flint_abort();
        }

        if (result)
        {
            for (i = 0; i < n; i++)
            {
                count = 0;
                for (j = 0; j < n; j++)
                    count += acb_contains(F + i, E + j);

                count2 = 0;
                for (j = 0; j < n; j++)
                    count2 += acb_overlaps(F + i, E + j);

                if (count != 1 || count2 != 1)
                {
                    flint_printf("FAIL: count\n\n");
                    flint_printf("algorithm = %d\n\n", algorithm);
                    flint_printf("A = \n"); acb_mat_printd(A, 20); flint_printf("\n\n");
                    flint_printf("R = \n"); acb_mat_printd(R, 20); flint_printf("\n\n");
                    flint_printf("L = \n"); acb_mat_printd(L, 20); flint_printf("\n\n");
                    flint_printf("D = \n"); acb_mat_printd(D, 20); flint_printf("\n\n");
                    flint_printf("LAR = \n"); acb_mat_printd(LAR, 20); flint_printf("\n\n");
                    flint_printf("i = %wd, count = %wd, count2 = %wd\n\n", i, count, count2);
                    flint_abort();
                }
            }
        }

        acb_mat_clear(A);
        acb_mat_clear(L);
        acb_mat_clear(R);
        acb_mat_clear(LAR);
        acb_mat_clear(D);
        acb_clear(b);
        _acb_vec_clear(E, n);
        _acb_vec_clear(F, n);
    }

    /* Test convergence, given companion matrices */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_mat_t A, R, QC;
        acb_ptr E;
        acb_ptr roots;
        fmpq_mat_t Q, Qinv;
        acb_poly_t f;
        slong i, j, n, prec, count, count2;
        int algorithm, success;

        algorithm = n_randint(state, 3);
        n = n_randint(state, 10);
        roots = _acb_vec_init(n);
        E = _acb_vec_init(n);
        acb_poly_init(f);
        acb_mat_init(A, n, n);
        acb_mat_init(R, n, n);
        fmpq_mat_init(Q, n, n);
        fmpq_mat_init(Qinv, n, n);
        acb_mat_init(QC, n, n);

        for (i = 0; i < n; i++)
        {
            new_root:
            acb_randtest(roots + i, state, 2 + n_randint(state, 100), 4);
            acb_get_mid(roots + i, roots + i);

            for (j = 0; j < i; j++)
                if (acb_equal(roots + i, roots + j))
                    goto new_root;
        }

        do {
            fmpq_mat_randtest(Q, state, 2 + n_randint(state, 100));
        } while (!fmpq_mat_inv(Qinv, Q));

        success = 0;

        for (prec = 32; !success; prec *= 2)
        {
            if (prec > 10000)
            {
                flint_printf("FAIL: unsuccessful, prec > 10000\n\n");
                flint_printf("algorithm = %d, iter %wd\n\n", algorithm, iter);
                flint_printf("A = \n"); acb_mat_printd(A, 20); flint_printf("\n\n");
                flint_printf("R = \n"); acb_mat_printd(R, 20); flint_printf("\n\n");
                flint_printf("roots = \n");
                for (j = 0; j < n; j++)
                {
                    acb_printd(roots + j, 20);
                    flint_printf("\n");
                }
                flint_abort();
            }

            acb_poly_product_roots(f, roots, n, prec);

            acb_mat_companion(A, f, prec);
            acb_mat_set_fmpq_mat(QC, Q, prec);
            acb_mat_mul(A, A, QC, prec);
            acb_mat_set_fmpq_mat(QC, Qinv, prec);
            acb_mat_mul(A, QC, A, prec);

            acb_mat_approx_eig_qr(E, NULL, R, A, NULL, 0, prec);

            if (algorithm == 0)
                success = acb_mat_eig_simple(E, NULL, NULL, A, E, R, prec);
            else if (algorithm == 1)
                success = acb_mat_eig_simple_rump(E, NULL, NULL, A, E, R, prec);
            else
                success = acb_mat_eig_simple_vdhoeven_mourrain(E, NULL, NULL, A, E, R, prec);

            if (success)
            {
                for (i = 0; i < n; i++)
                {
                    count = 0;
                    for (j = 0; j < n; j++)
                        count += acb_contains(E + i, roots + j);

                    count2 = 0;
                    for (j = 0; j < n; j++)
                        count2 += acb_overlaps(E + i, roots + j);

                    if (count != 1 || count2 != 1)
                    {
                        flint_printf("FAIL: count\n\n");
                        flint_printf("algorithm = %d\n\n", algorithm);
                        flint_printf("A = \n"); acb_mat_printd(A, 20); flint_printf("\n\n");
                        flint_printf("R = \n"); acb_mat_printd(R, 20); flint_printf("\n\n");
                        flint_printf("i = %wd, count = %wd, count2 = %wd\n\n", i, count, count2);
                        flint_printf("roots = \n");
                        for (j = 0; j < n; j++)
                        {
                            acb_printd(roots + j, 20);
                            flint_printf("\n");
                        }
                        flint_printf("E = \n");
                        for (j = 0; j < n; j++)
                        {
                            acb_printd(E + j, 20);
                            flint_printf("\n");
                        }
                        flint_abort();
                    }
                }
            }
        }

        fmpq_mat_clear(Q);
        fmpq_mat_clear(Qinv);
        acb_mat_clear(QC);
        acb_mat_clear(A);
        acb_mat_clear(R);
        acb_poly_clear(f);
        _acb_vec_clear(roots, n);
        _acb_vec_clear(E, n);
    }

    TEST_FUNCTION_END(state);
}
