/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int main(void)
{
    slong iter;
    flint_rand_t state;

    flint_printf("ql_reduce....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: agrees with naive algorithms */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 2 + n_randint(state, 3);
        slong n = 1 << g;
        slong prec = ACB_THETA_LOW_PREC + n_randint(state, 100);
        slong bits = 4;
        slong s = n_randint(state, g + 1);
        acb_mat_t tau, tau0;
        arb_mat_t Y;
        acb_ptr z, new_z, th, th0, test;
        arb_ptr x;
        acb_t c;
        arb_t u, abs;
        ulong a0, a1, b0, b1, fixed_a1;
        slong j, k;

        acb_mat_init(tau, g, g);
        arb_mat_init(Y, g, g);
        z = _acb_vec_init(g);
        new_z = _acb_vec_init(g);
        th = _acb_vec_init(n * n);
        th0 = _acb_vec_init(n * n);
        test = _acb_vec_init(n * n);
        x = _arb_vec_init(g);
        acb_init(c);
        arb_init(u);
        arb_init(abs);

        /* Make period matrix with splitting at s */
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        for (j = s; j < g; j++)
        {
            for (k = s; k < g; k++)
            {
                acb_mul_2exp_si(acb_mat_entry(tau, j, k),
                    acb_mat_entry(tau, j, k), 10);
            }
        }

        /* Choose z as Y.v + error with v either 0, 1/4 or 1/2 entries */
        for (k = 0; k < g; k++)
        {
            arb_set_si(&x[k], n_randint(state, 3));
        }
        _arb_vec_scalar_mul_2exp_si(x, x, g, -2);
        arb_mat_vector_mul_col(x, Y, x, prec);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
            arb_add(acb_imagref(&z[k]), acb_imagref(&z[k]), &x[k], prec);
        }

        s = acb_theta_ql_reduce(new_z, c, u, &fixed_a1, z, tau, prec);
        acb_theta_naive_all(th, z, 1, tau, prec);

        /* If s == -1, check that theta values are small */
        if (s == -1)
        {
            for (k = 0; k < n * n; k++)
            {
                acb_abs(abs, &th[k], prec);
                if (arb_gt(abs, u))
                {
                    flint_printf("FAIL (g = %wd, s = %wd)", g, s);
                    acb_mat_printd(tau, 5);
                    flint_printf("values, bound:\n");
                    _acb_vec_printd(th, n * n, 5);
                    arb_printd(u, 5);
                    flint_printf("\n");
                    flint_abort();
                }
            }
        }

        /* Otherwise, construct test vector */
        if (s == 0)
        {
            acb_one(&th0[0]);
        }
        else
        {
            acb_mat_window_init(tau0, tau, 0, s, 0, s);
            acb_theta_naive_all(th0, new_z, 1, tau0, prec);
            acb_mat_window_clear(tau0);
        }

        for (k = 0; k < n * n; k++)
        {
            a0 = k >> (g + g - s);
            a1 = (k >> g) % (1 << (g - s));
            b0 = k >> (g - s) % (1 << s);
            b1 = k % (1 << (g - s));
            if (a1 == fixed_a1)
            {
                acb_mul(&test[k], c, &th0[(a0 << s) + b0], prec);
                acb_mul_powi(&test[k], &test[k], acb_theta_char_dot(a1, b1, g));
            }
            acb_add_error_arb(&test[k], u);
        }

        if (!_acb_vec_overlaps(th, test, n * n))
        {
            flint_printf("FAIL (g = %wd, s = %wd)", g, s);
            acb_mat_printd(tau, 5);
            flint_printf("th, test:\n");
            _acb_vec_printd(th, n * n, 5);
            _acb_vec_printd(test, n * n, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        arb_mat_clear(Y);
        _acb_vec_clear(z, g);
        _acb_vec_clear(new_z, g);
        _acb_vec_clear(th, n * n);
        _acb_vec_clear(th0, n * n);
        _acb_vec_clear(test, n * n);
        _arb_vec_clear(x, g);
        acb_clear(c);
        arb_clear(u);
        arb_clear(abs);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}

