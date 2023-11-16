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

TEST_FUNCTION_START(acb_mat_eig_global_enclosure, state)
{
    slong iter;

    for (iter = 0; iter < 250 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_mat_t A, X;
        acb_ptr E, F, B;
        mag_t eps;
        slong i, j, n, prec;
        int success;

        n = n_randint(state, 8);
        prec = 2 + n_randint(state, 100);

        acb_mat_init(A, n, n);
        acb_mat_init(X, n, n);
        E = _acb_vec_init(n);
        F = _acb_vec_init(n);
        B = _acb_vec_init(n);
        mag_init(eps);

        for (i = 0; i < n; i++)
            acb_randtest(E + i, state, prec, 3);

        acb_mat_randtest_eig(A, state, E, prec);
        acb_mat_approx_eig_qr(F, NULL, X, A, NULL, 0, prec);

        /* perturb F further */
        if (n_randint(state, 2))
        {
            for (i = 0; i < n; i++)
            {
                acb_randtest(B, state, prec, 1);
                acb_mul_2exp_si(B, B, -n_randint(state, prec));
                acb_add(F + i, F + i, B, prec);
            }
        }

        acb_mat_eig_global_enclosure(eps, A, F, X, prec);

        _acb_vec_set(B, F, n);
        for (i = 0; i < n; i++)
            acb_add_error_mag(B + i, eps);

        for (i = 0; i < n; i++)
        {
            success = 0;

            for (j = 0; j < n; j++)
            {
                if (acb_contains(B + j, E + i))
                {
                    success = 1;
                    break;
                }
            }

            if (!success)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A =\n");
                acb_mat_printd(A, 20);
                flint_printf("\n\ni = %wd\n\n", i);
                flint_printf("\nE =\n");
                for (j = 0; j < n; j++)
                {
                    acb_printn(E + j, 20, 0); flint_printf("\n");
                }
                flint_printf("\nB =\n");
                for (j = 0; j < n; j++)
                {
                    acb_printn(B + j, 20, 0); flint_printf("\n");
                }
                flint_abort();
            }
        }

        acb_mat_clear(A);
        acb_mat_clear(X);
        _acb_vec_clear(E, n);
        _acb_vec_clear(F, n);
        _acb_vec_clear(B, n);
        mag_clear(eps);
    }

    TEST_FUNCTION_END(state);
}
