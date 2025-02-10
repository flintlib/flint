/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "arb.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_agm_mul_tight, state)
{
    slong iter;

    /* Test: respects relative precision, and agrees with agm_mul at huge precision */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 4);
        slong n = 1 << g;
        slong mprec = 50 + n_randint(state, 200);
        slong prec = mprec + 500;
        slong bits = n_randint(state, 3);
        slong delta = 25;
        int all = iter % 2;
        slong nbth = (all ? n * n : n);
        acb_mat_t tau;
        acb_ptr z;
        acb_ptr th, r, test;
        arb_ptr ds, exps;
        arb_t x, t;
        arf_t eps, u;
        slong k, b;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(2 * g);
        r = _acb_vec_init(nbth);
        test = _acb_vec_init(nbth);
        th = _acb_vec_init(2 * n);
        ds = _arb_vec_init(2 * n);
        exps = _arb_vec_init(2 * n);
        arb_init(x);
        arb_init(t);
        arf_init(eps);
        arf_init(u);

        /* Generate distances, not too crazy */
        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_randtest_vec_reduced(z + g, state, 1, tau, 0, prec);
        acb_theta_eld_distances(ds, z, 2, tau, prec);

        /* Generate values */
        for (k = 0; k < 2 * n; k++)
        {
            arb_neg(&exps[k], &ds[k]);
            arb_exp(&exps[k], &exps[k], prec);
            acb_urandom(&th[k], state, prec);
            acb_mul_arb(&th[k], &th[k], &exps[k], prec);
        }

        acb_theta_agm_mul(test, th, th + n, g, all, prec);

        /* Add error so that relative precision is roughly mprec */
        for (k = 0; k < 2 * n; k++)
        {
            arb_one(x);
            arb_mul_2exp_si(x, x, -mprec);
            arb_mul(x, x, &exps[k], prec);
            acb_add_arb(&th[k], &th[k], x, prec);
            acb_add_error_arb(&th[k], x);
        }

        acb_theta_agm_mul_tight(r, th, th + n, ds, ds + n, g, all, mprec);

        if (!_acb_vec_overlaps(r, test, nbth)
            || !_acb_vec_is_finite(r, nbth)
            || !_acb_vec_is_finite(test, nbth))
        {
            flint_printf("FAIL (overlap)\n");
            flint_abort();
        }

        for (k = 0; k < n; k++)
        {
            if (all)
            {
                arb_zero(x);
                for (b = 0; b < n; b++)
                {
                    acb_abs(t, &r[k * n + b], prec);
                    arb_max(x, x, t, prec);
                }
            }
            else
            {
                acb_abs(x, &r[k], prec);
            }
            arb_mul_2exp_si(x, x, -g);

            if (arb_gt(x, &exps[n + k]))
            {
                flint_printf("FAIL (absolute value, k = %wd)\n", k);
                flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
                acb_mat_printd(tau, 5);
                flint_printf("distances:\n");
                _arb_vec_printd(ds, 2 * n, 5);
                flint_printf("values:\n");
                _acb_vec_printd(th, 2 * n, 5);
                flint_printf("result:\n");
                _acb_vec_printd(r, nbth, 5);
                flint_abort();
            }

            if (all)
            {
                arf_zero(eps);
                for (b = 0; b < n; b++)
                {
                    acb_get_rad_ubound_arf(u, &r[k * n + b], prec);
                    arf_max(eps, eps, u);
                }
            }
            else
            {
                acb_get_rad_ubound_arf(eps, &r[k], prec);
            }
            arb_set_arf(x, eps);
            arb_mul_2exp_si(x, x, mprec - delta);

            if (arb_gt(x, &exps[n + k]))
            {
                flint_printf("FAIL (precision loss, k = %wd)\n", k);
                flint_printf("g = %wd, prec = %wd, tau:\n", g, prec);
                acb_mat_printd(tau, 5);
                flint_printf("distances:\n");
                _arb_vec_printd(ds, 2 * n, 5);
                flint_printf("values:\n");
                _acb_vec_printd(th, 2 * n, 5);
                flint_printf("result:\n");
                _acb_vec_printd(r, nbth, 5);
                flint_printf("x:\n");
                arb_printd(x, 5);
                flint_printf("\n");
                flint_abort();
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, 2 * g);
        _acb_vec_clear(r, nbth);
        _acb_vec_clear(test, nbth);
        _acb_vec_clear(th, 2 * n);
        _arb_vec_clear(ds, 2 * n);
        _arb_vec_clear(exps, 2 * n);
        arb_clear(x);
        arb_clear(t);
        arf_clear(eps);
        arf_clear(u);
    }

    TEST_FUNCTION_END(state);
}
