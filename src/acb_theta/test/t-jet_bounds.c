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

    flint_printf("jet_bounds....");
    fflush(stdout);

    flint_randinit(state);

    /* Test: bounds are finite, theta values correctly bounded */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong lp = ACB_THETA_LOW_PREC;
        slong prec = lp + n_randint(state, 1000);
        slong bits = n_randint(state, 4);
        slong ord = 1 + n_randint(state, 10);
        slong g = 1 + n_randint(state, 3);
        slong n2 = 1 << (2 * g);
        acb_mat_t tau;
        acb_ptr z, x, th;
        arb_t c, eps, rho, abs, t;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        x = _acb_vec_init(g);
        th = _acb_vec_init(n2);
        arb_init(c);
        arb_init(eps);
        arb_init(rho);
        arb_init(abs);
        arb_init(t);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_scalar_mul_2exp_si(tau, tau, -2);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }

        acb_theta_jet_bounds(eps, c, rho, z, tau, ord, prec, lp);
        arb_mul_2exp_si(c, c, prec - 1);

        if (!arb_is_finite(eps) || !arb_is_finite(c))
        {
            flint_printf("FAIL (infinite)\n");
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            flint_printf("eps, c, rho:\n");
            arb_printd(eps, 10);
            flint_printf("\n");
            arb_printd(c, 10);
            flint_printf("\n");
            arb_printd(rho, 10);
            flint_printf("\n");
            flint_abort();
        }

        for (k = 0; k < g; k++)
        {
            acb_urandom(&x[k], state, prec);
        }
        _acb_vec_scalar_mul_arb(x, x, g, rho, prec);
        _acb_vec_add(x, x, z, g, prec);
        acb_theta_naive_all(th, x, 1, tau, lp);

        arb_zero(abs);
        for (k = 0; k < n2; k++)
        {
            acb_abs(t, &th[k], lp);
            arb_max(abs, abs, t, lp);
        }

        if (arb_gt(abs, c))
        {
            flint_printf("FAIL (bound)\n");
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            flint_printf("rho, c, abs:\n");
            arb_printd(rho, 10);
            flint_printf("\n");
            arb_printd(c, 10);
            flint_printf("\n");
            arb_printd(abs, 10);
            flint_printf("\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(z, g);
        _acb_vec_clear(x, g);
        _acb_vec_clear(th, n2);
        arb_clear(c);
        arb_clear(eps);
        arb_clear(rho);
        arb_clear(abs);
        arb_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
