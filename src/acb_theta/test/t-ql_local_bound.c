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

TEST_FUNCTION_START(acb_theta_ql_local_bound, state)
{
    slong iter;

    /* Test: bounds are finite, theta values correctly bounded */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong lp = ACB_THETA_LOW_PREC;
        slong prec = lp + 100;
        slong bits = n_randint(state, 4);
        slong ord = 1 + n_randint(state, 10);
        slong n2 = 1 << (2 * g);
        acb_mat_t tau;
        acb_ptr z, x, th;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_t ctx;
        arb_ptr d;
        arb_t c, rho, abs, t;
        slong k;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        x = _acb_vec_init(g);
        th = _acb_vec_init(n2);
        acb_theta_ctx_tau_init(ctx_tau, 1, g);
        acb_theta_ctx_z_init(ctx, g);
        d = _arb_vec_init(n);
        arb_init(c);
        arb_init(rho);
        arb_init(abs);
        arb_init(t);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_mat_scalar_mul_2exp_si(tau, tau, -2);
        for (k = 0; k < g; k++)
        {
            acb_urandom(&z[k], state, prec);
        }

        acb_theta_ql_local_bound(c, rho, z, tau, ord);

        if (!arb_is_finite(rho)
            || !arb_is_finite(c)
            || arb_contains_zero(rho))
        {
            flint_printf("FAIL (infinite)\n");
            acb_mat_printd(tau, 5);
            _acb_vec_printd(z, g, 5);
            flint_printf("c, rho:\n");
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

        acb_theta_ctx_tau_set(ctx_tau, tau, lp);
        acb_theta_ctx_z_set(ctx, x, ctx_tau, lp);
        acb_theta_eld_distances(d, x, 1, tau, lp);
        acb_theta_sum(th, ctx, 1, ctx_tau, d, 1, 1, 1, lp);
        _acb_vec_scalar_mul_arb(th, th, n2, &ctx->u, lp);

        arb_zero(abs);
        for (k = 0; k < n2; k++)
        {
            acb_abs(t, &th[k], lp);
            arb_max(abs, abs, t, lp);
        }

        if (arb_gt(abs, c) || !arb_is_finite(abs))
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
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_clear(ctx);
        _arb_vec_clear(d, n);
        arb_clear(c);
        arb_clear(rho);
        arb_clear(abs);
        arb_clear(t);
    }

    TEST_FUNCTION_END(state);
}
