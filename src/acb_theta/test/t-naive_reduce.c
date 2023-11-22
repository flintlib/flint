/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_naive_reduce, state)
{
    slong iter;

    /* Test: special values of z */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 5);
        slong nbz = n_randint(state, 5);
        slong bits = n_randint(state, 5);
        slong prec = 100 + n_randint(state, 200);
        acb_mat_t tau;
        arb_mat_t Y, C;
        acb_ptr z, new_z, c;
        arb_ptr u, v, w, a;
        acb_t t, x;
        slong *n, *zero;
        slong err_exp = - 10 - n_randint(state, 20);
        slong k, j;
        int res;

        acb_mat_init(tau, g, g);
        arb_mat_init(Y, g, g);
        arb_mat_init(C, g, g);
        z = _acb_vec_init(g * nbz);
        new_z = _acb_vec_init(g * nbz);
        c = _acb_vec_init(nbz);
        u = _arb_vec_init(nbz);
        v = _arb_vec_init(g);
        w = _arb_vec_init(g * nbz);
        a = _arb_vec_init(g * nbz);
        acb_init(t);
        acb_init(x);
        n = flint_malloc(g * nbz * sizeof(slong));
        zero = flint_malloc(g * sizeof(slong));

        /* Set tau, cho, Y */
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_cho(C, tau, prec);
        acb_mat_get_imag(Y, tau);

        /* Test: if z is real, c = 1, u = 1 and v = 0 */
        for (k = 0; k < g * nbz; k++)
        {
            arb_randtest_precise(acb_realref(&z[k]), state, prec, bits);
        }
        acb_theta_naive_reduce(v, new_z, a, c, u, z, nbz, tau, prec);

        res = 1;
        for (k = 0; k < nbz; k++)
        {
            res = res && acb_is_one(&c[k]);
            res = res && arb_is_one(&u[k]);
        }

        if (!_arb_vec_is_zero(v, g) || !res)
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        /* Test: if im(z) = - Y . (even integral vector n) + small error,
           then terms for n and 0 correspond and v is small */
        for (j = 0; j < g; j++)
        {
            zero[j] = 0;
        }
        for (k = 0; k < nbz; k++)
        {
            for (j = k * g; j < (k + 1) * g; j++)
            {
                n[j] = 2 * n_randint(state, 10);
                arb_set_si(&w[j], n[j]);
            }
            arb_mat_vector_mul_col(w + k * g, Y, w + k * g, prec);
            for (j = k * g; j < (k + 1) * g; j++)
            {
                arb_urandom(acb_imagref(&z[j]), state, prec);
                arb_mul_2exp_si(acb_imagref(&z[j]), acb_imagref(&z[j]), err_exp);
                arb_sub(acb_imagref(&z[j]), acb_imagref(&z[j]), &w[j], prec);
            }
        }
        acb_theta_naive_reduce(v, new_z, a, c, u, z, nbz, tau, prec);

        for (k = 0; k < g * nbz; k++)
        {
            if (!arb_equal_si(&a[k], -n[k]))
            {
                flint_printf("FAIL (integral vector)\n");
                _arb_vec_printd(a, g * nbz, 5);
                flint_printf("k = %wd, n[k] = %wd\n", k, n[k]);
                flint_abort();
            }
        }

        for (k = 0; k < nbz; k++)
        {
            acb_theta_naive_term(x, z + k * g, tau, NULL, n + k * g, prec);
            acb_theta_naive_term(t, new_z + k * g, tau, NULL, zero, prec);
            acb_mul(t, t, &c[k], prec);

            if (!acb_overlaps(x, t))
            {
                flint_printf("FAIL (value, k = %wd)\n", k);
                flint_printf("tau:\n");
                acb_mat_printd(tau, 10);
                flint_printf("z:\n");
                _acb_vec_printd(z + k * g, g, 10);
                flint_printf("values:\n");
                acb_printd(x, 10);
                flint_printf("\n");
                acb_printd(t, 10);
                flint_printf("\n");
                acb_printd(&c[k], 10);
                flint_printf("\nNew z:\n");
                _acb_vec_printd(new_z + k * g, g, 10);
                flint_abort();
            }
        }

        arb_mat_inv(C, C, prec);
        arb_mat_vector_mul_col(v, C, v, prec);
        for (k = 0; k < g; k++)
        {
            arb_mul_2exp_si(&v[k], &v[k], - err_exp - 2);
            arb_sub_si(&v[k], &v[k], 1, prec);
            if (!arb_is_negative(&v[k]))
            {
                flint_printf("FAIL (offset)\n");
                arb_printd(&v[k], 10);
                flint_printf("\n");
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        arb_mat_clear(Y);
        arb_mat_clear(C);
        _acb_vec_clear(z, g * nbz);
        _acb_vec_clear(new_z, g * nbz);
        _acb_vec_clear(c, nbz);
        _arb_vec_clear(u, nbz);
        _arb_vec_clear(v, g);
        _arb_vec_clear(w, g * nbz);
        _arb_vec_clear(a, g * nbz);
        acb_clear(t);
        acb_clear(x);
        flint_free(n);
        flint_free(zero);
    }

    TEST_FUNCTION_END(state);
}

