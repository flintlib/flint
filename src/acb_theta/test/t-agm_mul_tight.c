/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_agm_mul_tight, state)
{
    slong iter;

    /* Test: respects relative precision */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong n = 1 << g;
        slong mprec = 50 + n_randint(state, 500);
        slong prec = mprec + 50;
        slong bits = n_randint(state, 3);
        slong delta = 25;
        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th, th0, r;
        arb_ptr ds;
        arb_t x, t;
        arf_t eps;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(2 * g);
        r = _acb_vec_init(n);
        th = _acb_vec_init(n);
        th0 = _acb_vec_init(n);
        ds = _arb_vec_init(2 * n);
        arb_init(x);
        arb_init(t);
        arf_init(eps);

        /* Generate distances, not too crazy */
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_randtest_vec_reduced(z + g, state, 1, tau, 0, prec);
        acb_theta_agm_distances(ds, z, 2, tau, prec);

        /* Generate values */
        for (k = 0; k < n; k++)
        {
            arb_neg(x, &ds[n + k]);
            arb_exp(x, x, prec);
            acb_urandom(&th[k], state, prec);
            acb_mul_arb(&th[k], &th[k], x, prec);

            arb_neg(x, &ds[k]);
            arb_exp(x, x, prec);
            acb_urandom(&th0[k], state, prec);
            acb_mul_arb(&th0[k], &th0[k], x, prec);
        }

        acb_theta_agm_mul_tight(r, th0, th, ds, ds + n, g, mprec);

        for (k = 0; k < n; k++)
        {
            acb_abs(x, &r[k], prec);
            arb_neg(t, &ds[n + k]);
            arb_exp(t, t, prec);
            if (arb_gt(x, t))
            {
                flint_printf("FAIL (absolute value, k = %wd)\n", k);
                flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
                acb_mat_printd(tau, 5);
                flint_printf("distances:\n");
                _arb_vec_printd(ds, n, 5);
                _arb_vec_printd(ds + n, n, 5);
                flint_printf("values:\n");
                _acb_vec_printd(th0, n, 5);
                _acb_vec_printd(th, n, 5);
                flint_printf("result:\n");
                _acb_vec_printd(r, n, 5);
                flint_abort();
            }

            acb_get_rad_ubound_arf(eps, &r[k], prec);
            arb_set_arf(x, eps);
            arb_mul_2exp_si(t, t, -mprec + delta);
            if (arb_gt(x, t))
            {
                flint_printf("FAIL (precision loss, k = %wd)\n", k);
                flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
                acb_mat_printd(tau, 5);
                flint_printf("distances:\n");
                _arb_vec_printd(ds, n, 5);
                _arb_vec_printd(ds + n, n, 5);
                flint_printf("values:\n");
                _acb_vec_printd(th0, n, 5);
                _acb_vec_printd(th, n, 5);
                flint_printf("result:\n");
                _acb_vec_printd(r, n, 5);
                flint_printf("x, t:\n");
                arb_printd(x, 5);
                flint_printf("\n");
                arb_printd(t, 5);
                flint_printf("\n");
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(r, n);
        _acb_vec_clear(th, n);
        _acb_vec_clear(th0, n);
        _arb_vec_clear(ds, 2 * n);
        arb_clear(x);
        arb_clear(t);
        arf_clear(eps);
    }

    TEST_FUNCTION_END(state);
}
