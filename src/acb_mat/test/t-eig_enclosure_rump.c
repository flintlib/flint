/*
    Copyright (C) 2018 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"

TEST_FUNCTION_START(acb_mat_eig_enclosure_rump, state)
{
    slong iter;

    /* Test random matrices */
    for (iter = 0; iter < 1000 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_mat_t A, X, R, AR, J, RJ, Z, Z0;
        acb_ptr E, F;
        acb_t b, lambda;
        slong i, j, h, k, n, iter2, prec, found_eigenvalue;

        n = 1 + n_randint(state, 7);
        prec = 2 + n_randint(state, 200);

        acb_mat_init(A, n, n);
        acb_mat_init(X, n, n);
        acb_init(lambda);
        acb_init(b);
        E = _acb_vec_init(n);
        F = _acb_vec_init(n);

        if (n_randint(state, 2))
        {
            for (i = 0; i < n; i++)
                acb_randtest(E + i, state, prec, 3);
        }
        else
        {
            /* Randomly repeat eigenvalues. */
            for (i = 0; i < n; i++)
            {
                if (i == 0 || n_randint(state, 2))
                    acb_randtest(E + i, state, prec, 3);
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
        acb_mat_approx_eig_qr(F, NULL, X, A, NULL, 0, prec);

        /* Perturb F further. */
        if (n_randint(state, 4) == 0)
        {
            for (i = 0; i < n; i++)
            {
                acb_randtest(b, state, prec, 1);
                acb_mul_2exp_si(b, b, -n_randint(state, prec));
                acb_add(F + i, F + i, b, prec);
            }
        }

        /* Perturb X further. */
        if (n_randint(state, 10) == 0)
        {
            j = n_randint(state, n);

            for (i = 0; i < n; i++)
            {
                acb_randtest(b, state, prec, 1);
                acb_mul_2exp_si(b, b, -10 - n_randint(state, prec));
                acb_add(acb_mat_entry(X, i, j), acb_mat_entry(X, i, j), b, prec);
            }
        }

        /* Test k = 1 */
        if (1)
        {
            acb_mat_init(R, n, 1);
            acb_mat_init(AR, n, 1);
            acb_mat_init(Z, n, 1);
            acb_mat_init(Z0, n, 1);

            for (j = 0; j < n; j++)
            {
                acb_set(lambda, F + j);
                for (i = 0; i < n; i++)
                    acb_set(acb_mat_entry(R, i, 0), acb_mat_entry(X, i, j));
                acb_mat_eig_enclosure_rump(lambda, NULL, R, A, lambda, R, prec);

                acb_mat_mul(AR, A, R, prec);
                acb_mat_neg(Z, AR);
                acb_mat_scalar_addmul_acb(Z, R, lambda, prec);

                if (!acb_mat_contains(Z, Z0))
                {
                    flint_printf("FAIL: not containing zero!\n\n");
                    flint_printf("A = \n"); acb_mat_printd(A, 20); flint_printf("\n\n");
                    flint_printf("R = \n"); acb_mat_printd(R, 20); flint_printf("\n\n");
                    flint_printf("lambda = \n"); acb_printd(lambda, 20); flint_printf("\n\n");
                    flint_printf("Z = \n"); acb_mat_printd(Z, 20); flint_printf("\n\n");
                    flint_printf("E = \n");
                    for (j = 0; j < n; j++)
                    {
                        acb_printd(E + j, 20);
                        flint_printf("\n");
                    }
                    flint_abort();
                }

                found_eigenvalue = 0;
                for (j = 0; j < n; j++)
                {
                    if (acb_contains(lambda, E + j))
                        found_eigenvalue++;
                }

                if (found_eigenvalue == 0)
                {
                    flint_printf("FAIL: eigenvalue not found\n\n");
                    flint_printf("A = \n"); acb_mat_printd(A, 20); flint_printf("\n\n");
                    flint_printf("R = \n"); acb_mat_printd(R, 20); flint_printf("\n\n");
                    flint_printf("lambda = \n"); acb_printd(lambda, 20); flint_printf("\n\n");
                    flint_printf("Z = \n"); acb_mat_printd(Z, 20); flint_printf("\n\n");
                    flint_printf("E = \n");
                    for (j = 0; j < n; j++)
                    {
                        acb_printd(E + j, 20);
                        flint_printf("\n");
                    }
                    flint_abort();
                }
            }

            acb_mat_clear(R);
            acb_mat_clear(AR);
            acb_mat_clear(Z);
            acb_mat_clear(Z0);
        }

        /* Test k > 1 */
        for (iter2 = 1; iter2 < n; iter2++)
        {
            k = n_randint(state, n + 1);
            k = FLINT_MAX(k, 2);

            acb_mat_init(R, n, k);
            acb_mat_init(AR, n, k);
            acb_mat_init(Z, n, k);
            acb_mat_init(Z0, n, k);
            acb_mat_init(J, k, k);
            acb_mat_init(RJ, n, k);

            /* Random selection */
            for (h = 0; h < k; h++)
            {
                j = n_randint(state, n);
                if (h == 0 || n_randint(state, 2))
                    acb_set(lambda, F + j);
                for (i = 0; i < n; i++)
                    acb_set(acb_mat_entry(R, i, h), acb_mat_entry(X, i, j));
            }

            acb_mat_eig_enclosure_rump(lambda, J, R, A, lambda, R, prec);

            /* AY = YJ */
            acb_mat_mul(AR, A, R, prec);
            acb_mat_mul(RJ, R, J, prec);
            acb_mat_sub(Z, AR, RJ, prec);

            if (!acb_mat_contains(Z, Z0))
            {
                flint_printf("FAIL: not containing zero! (k = %wd, prec = %wd)\n\n", k, prec);
                flint_printf("A = \n"); acb_mat_printd(A, 20); flint_printf("\n\n");
                flint_printf("R = \n"); acb_mat_printd(R, 20); flint_printf("\n\n");
                flint_printf("lambda = \n"); acb_printd(lambda, 20); flint_printf("\n\n");
                flint_printf("J = \n"); acb_mat_printd(J, 20); flint_printf("\n\n");
                flint_printf("Z = \n"); acb_mat_printd(Z, 20); flint_printf("\n\n");
                flint_printf("E = \n");
                for (j = 0; j < n; j++)
                {
                    acb_printd(E + j, 20);
                    flint_printf("\n");
                }
                flint_abort();
            }

            found_eigenvalue = 0;
            for (j = 0; j < n; j++)
            {
                if (acb_contains(lambda, E + j))
                    found_eigenvalue++;
            }

            if (found_eigenvalue < k)
            {
                flint_printf("FAIL: eigenvalue not found (k = %wd, found = %wd)\n\n", k, found_eigenvalue);
                flint_printf("A = \n"); acb_mat_printd(A, 20); flint_printf("\n\n");
                flint_printf("R = \n"); acb_mat_printd(R, 20); flint_printf("\n\n");
                flint_printf("lambda = \n"); acb_printd(lambda, 20); flint_printf("\n\n");
                flint_printf("Z = \n"); acb_mat_printd(Z, 20); flint_printf("\n\n");
                flint_printf("E = \n");
                for (j = 0; j < n; j++)
                {
                    acb_printd(E + j, 20);
                    flint_printf("\n");
                }
                flint_abort();
            }

            acb_mat_clear(R);
            acb_mat_clear(AR);
            acb_mat_clear(Z);
            acb_mat_clear(Z0);
            acb_mat_clear(J);
            acb_mat_clear(RJ);
        }

        acb_mat_clear(A);
        acb_mat_clear(X);
        acb_clear(lambda);
        acb_clear(b);
        _acb_vec_clear(E, n);
        _acb_vec_clear(F, n);
    }

    TEST_FUNCTION_END(state);
}
