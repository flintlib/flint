/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_ql_a0, state)
{
    slong iter;

    /* Test: agrees with ql_a0_naive */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong prec = (g > 1 ? 200 : 500) + n_randint(state, 500);
        slong bits = n_randint(state, 5);
        slong hprec = prec + 100;
        int hast = iter % 2;
        int hasz = (iter % 4) / 2;
        slong nbt = (hast ? 3 : 1);
        slong nbz = (hasz ? 2 : 1);
        slong guard = ACB_THETA_LOW_PREC;
        slong lp = ACB_THETA_LOW_PREC;
        acb_mat_t tau;
        acb_ptr z, t, r, test;
        arb_ptr d, d0;
        arb_t y;
        slong k;
        int res;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        t = _acb_vec_init(g);
        r = _acb_vec_init(nbz * nbt * n);
        test = _acb_vec_init(nbz * nbt * n);
        d = _arb_vec_init(n);
        d0 = _arb_vec_init(n);
        arb_init(y);

        res = 0;
        while(!res)
        {
            acb_siegel_randtest_reduced(tau, state, prec, bits);
            arb_sub_si(y, acb_imagref(acb_mat_entry(tau, g - 1, g - 1)), 200, prec);
            res = arb_is_negative(y);
        }

        if (hast)
        {
            for (k = 0; k < g; k++)
            {
                arb_urandom(acb_realref(&t[k]), state, hprec);
            }
        }

        acb_theta_dist_a0(d0, z, tau, lp);
        if (hasz)
        {
            acb_siegel_randtest_vec(z, state, g, prec);
        }
        acb_theta_dist_a0(d, z, tau, lp);

        res = acb_theta_ql_a0(r, t, z, d0, d, tau, guard, prec);
        acb_theta_ql_a0_naive(test, t, z, d0, d, tau, guard, hprec);

        if (res && !_acb_vec_overlaps(r, test, nbz * nbt * n))
        {
            flint_printf("FAIL\n");
            flint_printf("g = %wd, prec = %wd, hprec = %wd, hasz = %wd, hast = %wd, tau:\n",
                g, prec, hprec, hasz, hast);
            acb_mat_printd(tau, 5);
            flint_printf("output:\n");
            _acb_vec_printd(r, nbz * nbt * n, 5);
            _acb_vec_printd(test, nbz * nbt * n, 5);
            flint_printf("difference:\n");
            _acb_vec_sub(r, r, test, nbz * nbt * n, hprec);
            _acb_vec_printd(r, nbz * nbt * n, 5);
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(t, g);
        _acb_vec_clear(r, nbz * nbt * n);
        _acb_vec_clear(test, nbz * nbt * n);
        _arb_vec_clear(d, n);
        _arb_vec_clear(d0, n);
        arb_clear(y);
    }

    TEST_FUNCTION_END(state);
}
