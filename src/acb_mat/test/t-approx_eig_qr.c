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

TEST_FUNCTION_START(acb_mat_approx_eig_qr, state)
{
    slong iter;

    /* Test random & DFT matrices */
    for (iter = 0; iter < 200 * 0.1 * flint_test_multiplier(); iter++)
    {
        acb_mat_t A, L, R;
        acb_ptr E;
        acb_t t;
        mag_t b;
        slong i, j, n, prec, goal, c0, c1, c2, c3;
        int wantL, wantR, result, dft;

        dft = n_randint(state, 2);
        if (dft)
            n = n_randint(state, 30);
        else
            n = n_randint(state, 15);
        goal = 2 + n_randint(state, 100);
        wantL = n_randint(state, 2);
        wantR = n_randint(state, 2);

        acb_mat_init(A, n, n);
        acb_mat_init(L, n, n);
        acb_mat_init(R, n, n);
        acb_init(t);
        mag_init(b);
        E = _acb_vec_init(n);

        for (prec = 32; ; prec *= 2)
        {
            if (dft)
            {
                acb_mat_dft(A, 0, prec);
            }
            else
            {
                acb_mat_randtest(A, state, 2 + n_randint(state, 200), 5);
                acb_mat_get_mid(A, A);
            }

            acb_mat_approx_eig_qr(E, wantL ? L : NULL, wantR ? R : NULL, A, NULL, 0, prec);

            if (dft)
            {
                /* Verify the known eigenvalues + multiplicities */
                c0 = c1 = c2 = c3 = 0;

                for (i = 0; i < n; i++)
                {
                    acb_set_d_d(t, 1.0, 0.0);
                    acb_sub(t, t, E + i, prec);
                    acb_get_mag(b, t);
                    c0 += (mag_cmp_2exp_si(b, -goal) < 0);

                    acb_set_d_d(t, -1.0, 0.0);
                    acb_sub(t, t, E + i, prec);
                    acb_get_mag(b, t);
                    c1 += (mag_cmp_2exp_si(b, -goal) < 0);

                    acb_set_d_d(t, 0.0, 1.0);
                    acb_sub(t, t, E + i, prec);
                    acb_get_mag(b, t);
                    c2 += (mag_cmp_2exp_si(b, -goal) < 0);

                    acb_set_d_d(t, 0.0, -1.0);
                    acb_sub(t, t, E + i, prec);
                    acb_get_mag(b, t);
                    c3 += (mag_cmp_2exp_si(b, -goal) < 0);
                }

                result = (n == 0 || (c0 == (n+4)/4 && c1 == (n+2)/4 && c2 == (n-1)/4 && c3 == (n+1)/4));
            }
            else
            {
                result = 1;
            }

            if (result && wantL)
            {
                acb_mat_t LA, D;
                acb_mat_init(LA, n, n);
                acb_mat_init(D, n, n);

                /* Check LA - lambda L = 0 */
                acb_mat_approx_mul(LA, L, A, prec);
                for (i = 0; i < n; i++)
                    acb_set(acb_mat_entry(D, i, i), E + i);
                acb_mat_approx_mul(D, D, L, prec);
                acb_mat_sub(LA, LA, D, prec);

                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        acb_get_mag(b, acb_mat_entry(LA, i, j));
                        result = result && (mag_cmp_2exp_si(b, -goal) < 0);
                    }
                }

                acb_mat_clear(LA);
                acb_mat_clear(D);
            }

            if (result && wantR)
            {
                acb_mat_t AR, D;
                acb_mat_init(AR, n, n);
                acb_mat_init(D, n, n);

                /* Check AR - R lambda = 0 */
                acb_mat_approx_mul(AR, A, R, prec);
                for (i = 0; i < n; i++)
                    acb_set(acb_mat_entry(D, i, i), E + i);
                acb_mat_approx_mul(D, R, D, prec);
                acb_mat_sub(AR, AR, D, prec);

                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        acb_get_mag(b, acb_mat_entry(AR, i, j));
                        result = result && (mag_cmp_2exp_si(b, -goal) < 0);
                    }
                }

                acb_mat_clear(AR);
                acb_mat_clear(D);
            }

            if (result)
                break;

            if (prec > 2000)
            {
                flint_printf("FAIL (convergence, dft = %d)\n\n", dft);
                flint_printf("n = %wd\n\n", n);
                acb_mat_printd(A, 10);
                flint_printf("\n\n");
                for (i = 0; i < n; i++)
                {
                    acb_printn(E + i, 50, 0);
                    flint_printf("\n");
                }
                flint_printf("\n");
                if (wantL)
                {
                    flint_printf("L = \n");
                    acb_mat_printd(L, 10);
                    flint_printf("\n\n");
                }
                if (wantR)
                {
                    flint_printf("R = \n");
                    acb_mat_printd(R, 10);
                    flint_printf("\n\n");
                }
                flint_abort();
            }
        }

        acb_mat_clear(A);
        acb_mat_clear(L);
        acb_mat_clear(R);
        _acb_vec_clear(E, n);
        acb_clear(t);
        mag_clear(b);
    }

    TEST_FUNCTION_END(state);
}
